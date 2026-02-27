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

# Точне знаходження ядра матриці через RREF
function exact_nullspace(A::Matrix{T}) where T <: Rational
    # Приводимо матрицю до ступенястого вигляду (RREF)
    R = RowEchelon.rref(A)
    
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

# Оптимізований solve_centralizer
function solve_centralizer(D::Derivation, max_deg::Int)
    vars = D.vars
    n = length(vars)
    
    # 1. Створюємо невідомі коефіцієнти
    # Використовуємо вектор для швидкого доступу
    coeffs = Num[]
    u_polys = [Num(0) for _ in 1:n]
    
    for p_idx in 1:n
        for d in 0:max_deg
            for i in 0:d
                j = d - i
                c_name = Symbol("c_$(p_idx)_$(i)_$(j)")
                c = first(@variables $c_name)
                push!(coeffs, c)
                u_polys[p_idx] += c * vars[1]^i * vars[2]^j
            end
        end
    end
    
    # 2. Обчислюємо комутатор
    Du = Derivation(u_polys, vars)
    bracket = lie_bracket(D, Du)
    
    # 3. Швидке формування системи рівнянь
    all_eqs = Num[]
    for p in bracket.polys
        # polynomial_coeffs вже повертає коефіцієнти при x^a y^b
        # Оскільки Du лінійний за coeffs, ці коефіцієнти будуть лінійними формами від c_...
        coeffs_dict, _ = Symbolics.polynomial_coeffs(p, vars)
        append!(all_eqs, values(coeffs_dict))
    end

    # 4. Будуємо матрицю A
    # Використання sparse допоможе на великих степенях
    A_sym = Symbolics.jacobian(all_eqs, coeffs)

    # 5. Створюємо розріджену числову матрицю
    # Використовуємо SparseMatrixCSC для економії пам'яті
    raw_vals = Symbolics.value.(A_sym)
    A_sparse = sparse(Matrix{Rational{Int128}}(raw_vals))
    
    # --- Блок аналізу розрідженості для PhD-звіту ---
    rows, cols = size(A_sparse)
    total_elements = rows * cols
    non_zeros = nnz(A_sparse)
    # Розрахунок відсотка нулів (sparsity)
    sparsity_percent = (1 - non_zeros / total_elements) * 100
    
    println("\n--- Статистика системи (max_deg = $max_deg) ---")
    println("Розмірність матриці: $rows x $cols")
    println("Ненульових елементів (nnz): $non_zeros")
    println("Степінь розрідженості (відсоток нулів): $(round(sparsity_percent, digits=2))%")
    # -----------------------------------------------
    
    # Перетворюємо на числову матрицю
    # Використовуємо BigInt для точності
    raw_vals = Symbolics.value.(A_sym)
    A_numeric = Matrix{Rational{Int128}}(raw_vals)
    # A_numeric = sparse(Matrix{Rational{BigInt}}(Symbolics.value.(A_sym)))
    # ns = nullspace(A_numeric)
    
    # 5. Знаходимо ядро (використовуємо вашу перевірену функцію)
    ns = exact_nullspace(A_numeric)
    
    return ns, coeffs, u_polys
end

end # module