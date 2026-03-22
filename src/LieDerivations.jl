module LieDerivations

using Symbolics
using LinearAlgebra
using SparseArrays
using RowEchelon


export Derivation, apply_derivation, lie_bracket, solve_centralizer

# Структура для диференціювання в W_n(K)
struct Derivation
    polys::Vector{Num}    # Поліноми P_i (коефіцієнти при ∂/∂x_i)
    vars::Vector{Num}     # Змінні x_1, ..., x_n
end

# Застосування D до функції f: D(f) = Σ P_i * ∂f/∂x_i
function apply_derivation(D::Derivation, f::Num)
    # Створюємо оператори диференціювання
    ops = [Differential(v)(f) for v in D.vars]
    # Обчислюємо результат символьно
    return expand_derivatives(sum(D.polys .* ops))
end

# Дужка Лі: оптимізована версія без проміжних спрощень
function lie_bracket(D1::Derivation, D2::Derivation)
    n = length(D1.vars)
    new_polys = Vector{Num}(undef, n)
    for i in 1:n
        # Обчислюємо компоненту комутатора
        # Прибираємо simplify всередині циклу
        res = apply_derivation(D1, D2.polys[i]) - apply_derivation(D2, D1.polys[i])
        new_polys[i] = expand(res) 
    end
    return Derivation(new_polys, D1.vars)
end

# Допоміжна функція для перевірки пропорційності двох диференціювань
function is_proportional(D1::Derivation, D2::Derivation)
    # Перевіряємо відношення компонентів: P1_i / P2_i має бути константою λ
    # Для PhD-точності краще перевіряти умову: P1_i * P2_j - P1_j * P2_i == 0
    vars = D1.vars
    n = length(vars)
    
    # Співвідношення між компонентами: P1[i]*P2[j] == P1[j]*P2[i]
    # Це еквівалентно пропорційності векторів у кожній точці
    for i in 1:n
        for j in (i+1):n
            cross_prod = simplify(D1.polys[i] * D2.polys[j] - D1.polys[j] * D2.polys[i])
            if !iszero(cross_prod)
                return false # Не пропорційні
            end
        end
    end
    return true # Пропорційні
end

# Універсальна функція для знаходження ядра
function exact_nullspace(A::AbstractMatrix{T}) where T <: Rational
    # Якщо матриця розріджена, конвертуємо її в Matrix лише для RowEchelon,
    # оскільки пакет RowEchelon працює з щільними масивами.
    # Проте ми зберігаємо інтерфейс для майбутніх розріджених оптимізацій.
    R = RowEchelon.rref(Matrix(A))
    
    m, n = size(R)
    pivot_cols = Int[]
    for i in 1:m
        j = findfirst(!iszero, R[i, :])
        if j !== nothing
            push!(pivot_cols, j)
        end
    end
    
    free_cols = setdiff(1:n, pivot_cols)
    null_basis = zeros(T, n, length(free_cols))
    
    for (k, j) in enumerate(free_cols)
        null_basis[j, k] = one(T)
        for i in 1:length(pivot_cols)
            null_basis[pivot_cols[i], k] = -R[i, j]
        end
    end
    return null_basis
end

function exact_sparse_nullspace(A::SparseMatrixCSC{T}) where T <: Rational
    m, n = size(A)
    S = copy(transpose(A)) 
    dropzeros!(S)

    pivot_rows = Int[] 
    pivot_cols = Int[] 
    
    # 1. Прямий хід
    for j in 1:m 
        nz_indices = S.colptr[j]:(S.colptr[j+1]-1)
        if isempty(nz_indices) continue end
        
        row_indices = S.rowval[nz_indices]
        available_rows = filter(r -> !(r ∈ pivot_cols), row_indices)
        
        if !isempty(available_rows)
            # Правило Марковіца
            i = available_rows[argmin([nnz(S[r, :]) for r in available_rows])]
            
            push!(pivot_cols, i)
            push!(pivot_rows, j)
            
            pivot_val = S[i, j]
            
            rows_to_update = findall(!iszero, S[i, :])
            for r in setdiff(rows_to_update, [j])
                factor = S[i, r] / pivot_val
                S[:, r] -= factor * S[:, j]
            end
            dropzeros!(S) 
        end
    end

    # --- ВИПРАВЛЕННЯ: Визначаємо змінні ПЕРЕД використанням ---
    free_vars = setdiff(1:n, pivot_cols)
    num_free = length(free_vars)
    
    # Створюємо порожню розріджену матрицю для базису
    null_basis = spzeros(T, n, num_free)
    
    # println("-->2 (Знайдено вільних змінних: $num_free)")

    # 2. Побудова базису ядра
    for (k, f_var) in enumerate(free_vars)
        # Встановлюємо 1 для вільної змінної
        null_basis[f_var, k] = one(T)
        # Обчислюємо значення для опорних змінних
        for (idx, p_var) in enumerate(pivot_cols)
            p_row = pivot_rows[idx]
            # Використовуємо S для знаходження коефіцієнтів
            if !iszero(S[f_var, p_row])
                null_basis[p_var, k] = -S[f_var, p_row] / S[p_var, p_row]
            end
        end
    end
    
    
    # println("--> Формування фінальної матриці")
    
    # 1. Витягуємо ненульові елементи (I - рядки, J - стовпці, V - значення)
    I, J, V = findnz(null_basis)
    
    # 2. Фільтруємо справжні нулі вручну (це надійніше за dropzeros! у складних випадках)
    # Залишаємо тільки ті індекси, де значення Rational не є нулем
    mask = .!iszero.(V)
    I_clean = I[mask]
    J_clean = J[mask]
    V_clean = V[mask]
    
    # 3. Створюємо абсолютно нову SparseMatrixCSC з чистими типами
    # Вказуємо тип T явно (це Rational{Int128})
    final_ns = sparse(I_clean, J_clean, V_clean, size(null_basis, 1), size(null_basis, 2))
    
    # println("--> Базис сформовано. Кількість елементів: $(length(V_clean))")
    
    return final_ns

end
    


function solve_centralizer(D::Derivation, max_deg::Int)
    vars = D.vars
    n = length(vars)
    k = 0 
    # Початковий стан результату
    # best_res = (ns=nothing, coeffs=nothing, u_polys=nothing, 
    #             is_proportional=false, is_valid=false, degree=0)
    sparsity_threshold = 90.0 # Поріг у відсотках

    all_found = []
    
    while k <= max_deg
        println("\n>>> Перевірка степеня k = $k")
        
        # 1. Генерація невідомих
        coeffs = Num[]
        u_polys = [Num(0) for _ in 1:n]
        for p_idx in 1:n
            for d in 0:k
                for i in 0:d
                    j = d - i
                    c_name = Symbol("c_$(p_idx)_$(i)_$(j)")
                    c = first(@variables $c_name)
                    push!(coeffs, c)
                    u_polys[p_idx] += c * vars[1]^i * vars[2]^j
                end
            end
        end
        
        # 2. Формування системи
        Du_sym = Derivation(u_polys, vars)
        bracket = lie_bracket(D, Du_sym)
        
        all_eqs = Num[]
        for p in bracket.polys
            coeffs_dict, _ = Symbolics.polynomial_coeffs(p, vars)
            append!(all_eqs, values(coeffs_dict))
        end

        # 4. Побудова символьної матриці Якобі
        A_sym = Symbolics.jacobian(all_eqs, coeffs)

        # 5. Глибоке очищення типів (PhD-Grade Cleaning)
        # Створюємо чистий масив раціональних чисел, ігноруючи символьні обгортки
        rows_dim, cols_dim = size(A_sym)
        A_numeric_clean = Matrix{Rational{Int128}}(undef, rows_dim, cols_dim)
        
        for c in 1:cols_dim
            for r in 1:rows_dim
                # Витягуємо значення, розгортаємо його і примусово кастуємо до Rational
                val = Symbolics.unwrap(Symbolics.value(A_sym[r, c]))
                # Якщо раптом потрапив Num(0) або подібне, беремо тільки числову частину
                A_numeric_clean[r, c] = Rational{Int128}(val)
            end
        end

        # Тепер створюємо розріджену матрицю з чистих чисел
        A_sparse = sparse(A_numeric_clean)
        dropzeros!(A_sparse)
        
        # Аналіз розрідженості
        non_zeros = nnz(A_sparse)
        sparsity_percent = (1 - non_zeros / (rows_dim * cols_dim)) * 100
        println("Розмірність: $rows_dim x $cols_dim, nnz: $non_zeros")
        println("Розрідженість: $(round(sparsity_percent, digits=2))%")
        
        # 6. Адаптивне знаходження ядра
        local ns
        if sparsity_percent > sparsity_threshold
            println("Метод: СПРАВЖНІЙ РОЗРІДЖЕНИЙ")
            ns = exact_sparse_nullspace(A_sparse)
        else
            println("Метод: КЛАСИЧНИЙ (DENSE)")
            ns = exact_nullspace(A_numeric_clean)
        end
        # println("result obtained, analyzing...")
        # 4. Аналіз результатів (ядро)
        if size(ns, 2) > 0
            for col in 1:size(ns, 2)
                # Конвертація у щільний вектор
                dense_col_vals = Vector{Rational{Int128}}(ns[:, col])
                sol_dict = Dict(coeffs[i] => dense_col_vals[i] for i in 1:length(coeffs))
                
                final_polys = [substitute(p, sol_dict) for p in u_polys]
                

                # 1. Пропусткаємо нульові розв'язки
                if all(iszero, final_polys)
                    continue
                end
                
                
                
                # 2. ПЕРЕВІРКА НА ВАЛІДНІСТЬ (Комутатор [D, Du] == 0)
                # Створюємо об'єкт диференціювання для знайденого розв'язку
                Du_found = Derivation(final_polys, vars)
                
                # Обчислюємо комутатор. Результуючі поліноми мають бути нульовими.
                comm_polys = lie_bracket(D, Du_found)
                is_valid = all(p -> iszero(simplify(p)), comm_polys.polys)

                if !is_valid
                    @warn "Знайдено розв'язок на k=$k, який не обнуляє комутатор! Пропускаємо."
                    continue
                end

                proportional = is_proportional(D, Du_found)

                is_interesting = !proportional

                # Додаємо у список тільки валідні та нетривіальні розв'язки
                push!(all_found, (
                    degree = k, 
                    polys = final_polys, 
                    is_interesting = !proportional,
                    is_valid = true
                ))
                
                
                status_type = !proportional ? "НЕ ПРОПОРЦІЙНИЙ" : "пропорційний"
                println("--> Збережено валідний $status_type розв'язок степеня $k")
            end
        end
        
        # Перехід до наступного степеня без виходу з функції
        k += 1
    end # кінець циклу while

    # 5. Повертаємо останній збережений "найкращий" результат
    return all_found
end

end # module