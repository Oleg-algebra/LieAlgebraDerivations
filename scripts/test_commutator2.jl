using LieDerivations
using Symbolics
using Latexify
using LinearAlgebra
using Revise
using LaTeXStrings 
using Dates       

# 1. Ініціалізація змінних
@variables x y

# 2. Визначення диференціювання
poly_x = y
poly_y = x
D = Derivation([poly_x, poly_y], [x, y])

max_deg = 15
println("Шукаємо централізатор для D (покроковий пошук до степеня $max_deg)...")

# "Розігрів" для JIT-компіляції (на малому степені)
_ = solve_centralizer(D, 1)

# 3. Основний виклик з вимірюванням часу
res = nothing
t_exec = @elapsed begin
    global res = solve_centralizer(D, 7)
end

println("Час виконання: $(round(t_exec, digits=3)) сек.")

println("Всього знайдено розв'язків: ", length(res))

# Вивести тільки цікаві (непропорційні)
for item in res
    if item.is_interesting
        println("Степінь $(item.degree): $(item.polys)")
    end
end



# Параметри вибору: "min" або "max"
res_type = "min" 

# 4. Форматування та запис результатів
filename = "tex-code/results.tex"

# 1. Функція для обчислення сумарного степеня поліномів (P1 + P2)
function get_total_poly_degree(polys)
    return sum(p -> iszero(p) ? 0 : Symbolics.degree(p), polys)
end

# 2. Попередня класифікація результатів та збір статистики
# all_found — це масив кортежів, який повернув ваш солвер
processed_results = []
stats_table = Dict{Int, Int}() # k => кількість цікавих (непропорційних)

for item in res
    if !item.is_valid continue end
    
    # Оновлюємо статистику для таблиці (тільки непропорційні)
    k = item.degree
    stats_table[k] = get(stats_table, k, 0) + (item.is_interesting ? 1 : 0)
    
    push!(processed_results, (
        polys = item.polys,
        k_degree = k,
        is_prop = !item.is_interesting, # цікавий = непропорційний
        sum_deg = get_total_poly_degree(item.polys)
    ))
end

# 3. Відбір згідно з параметром res_type (min або max)
# Пріоритет: спочатку шукаємо серед непропорційних
interesting = filter(r -> !r.is_prop, processed_results)
candidates = isempty(interesting) ? filter(r -> r.is_prop, processed_results) : interesting

final_selection = []
if !isempty(candidates)
    all_sum_degs = [c.sum_deg for c in candidates]
    target_val = (res_type == "min") ? minimum(all_sum_degs) : maximum(all_sum_degs)
    
    # Вибираємо всіх, хто має такий самий екстремальний сумарний степінь
    final_selection = filter(c -> c.sum_deg == target_val, candidates)
end

# 4. Запис у LaTeX файл
open(filename, "w") do io
    write(io, "% !TeX root = main.tex\n")
    write(io, "\\section*{Результати аналізу централізатора}\n")
    write(io, "\\noindent\\textit{Критерій вибору: $res_type сумарний степінь поліномів.}\\\\\n")

    # --- ТАБЛИЦЯ РОЗПОДІЛУ РОЗМІРНОСТЕЙ ---
    write(io, "\\subsection*{Розподіл нетривіальних елементів}\n")
    if !isempty(stats_table)
        sorted_k = sort(collect(keys(stats_table)))
        write(io, "\\begin{table}[h!]\n\\centering\n\\begin{tabular}{|l|" * "c|"^length(sorted_k) * "}\n\\hline\n")
        write(io, "Степінь \$ k \$ & " * join(sorted_k, " & ") * " \\\\ \\hline\n")
        write(io, "К-сть \$ D_u \\notin \\mathbb{K}D \$ & " * join([stats_table[k] for k in sorted_k], " & ") * " \\\\ \\hline\n")
        write(io, "\\end{tabular}\n\\caption{Розподіл розмірностей нетривіальної частини централізатора}\n\\end{table}\n\n")
    end

    # --- СПИСОК ОБРАНИХ ЕЛЕМЕНТІВ ---
    if isempty(final_selection)
        write(io, "Не знайдено жодного нетривіального розв'язку в заданому просторі.\n")
    else
        status_label = final_selection[1].is_prop ? "пропорційний" : "НЕ ПРОПОРЦІЙНИЙ"
        write(io, "\\subsection*{Обрані розв'язки ($status_label)}\n")
        write(io, "Ці елементи мають екстремальний сумарний степінь \$\\sum \\deg(P_i) = $(final_selection[1].sum_deg)\$.\\\\\n")
        
        write(io, "\\begin{itemize}\n")
        for c in final_selection
            # Факторизація для кращого вигляду в LaTeX
            factored_x = Symbolics.factorize(simplify(c.polys[1]))
            factored_y = Symbolics.factorize(simplify(c.polys[2]))
            
            l_x = latexify(factored_x, env=:raw)
            l_y = latexify(factored_y, env=:raw)
            
            write(io, "  \\item \$ D_{u} = ( $l_x ) \\partial_x + ( $l_y ) \\partial_y \$ \\quad (знайдено на \$ k = $(c.k_degree) \$)\n")
            
            # Вивід у консоль для швидкого контролю
            println("[$res_type] k=$(c.k_degree), sum_deg=$(c.sum_deg): D_u знайдено.")
        end
        write(io, "\\end{itemize}\n")
    end
    write(io, "\\vfill\\hrule\n")
end

println("\nЗвіт успішно згенеровано: $filename")