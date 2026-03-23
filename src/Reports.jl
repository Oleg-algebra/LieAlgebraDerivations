module Reports

using Symbolics
using Latexify
using Dates
using LinearAlgebra

export generate_latex_report

"""
Обчислює сумарний степінь координатних поліномів диференціювання.
"""
function get_total_poly_degree(polys)
    return sum(p -> iszero(p) ? 0 : Symbolics.degree(p), polys)
end

"""
Генерує комплексний LaTeX звіт.
- results: список all_found від солвера
- D_orig: початкове диференціювання
- vars: змінні [x, y]
- filename: шлях до .tex файлу
- res_type: "min" або "max" (критерій за сумарним степенем)
"""
function generate_latex_report(results, D_orig, vars, filename::String; res_type::String="min")
    # 1. Попередня обробка та статистика
    processed = []
    stats_table = Dict{Int, Int}() # k => к-сть непропорційних
    valid_count = 0
    invalid_count = 0

    for item in results
        if !item.is_valid
            invalid_count += 1
            continue 
        end
        valid_count += 1
        
        k = item.degree
        # Оновлюємо статистику для таблиці
        stats_table[k] = get(stats_table, k, 0) + (item.is_interesting ? 1 : 0)
        
        push!(processed, (
            polys = item.polys,
            k = k,
            is_prop = !item.is_interesting,
            sum_deg = get_total_poly_degree(item.polys)
        ))
    end

    # 2. Відбір за типом (min/max)
    # Пріоритет: непропорційні -> пропорційні
    interesting = filter(r -> !r.is_prop, processed)
    candidates = isempty(interesting) ? filter(r -> r.is_prop, processed) : interesting
    
    final_selection = []
    if !isempty(candidates)
        all_degs = [c.sum_deg for c in candidates]
        target = (res_type == "min") ? minimum(all_degs) : maximum(all_degs)
        final_selection = filter(c -> c.sum_deg == target, candidates)
    end

    # 3. Запис у LaTeX
    open(filename, "a") do io
        write(io, "%% Звіт згенеровано модулем Reports.jl\n")
        write(io, "%% Дата: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n\n")
        
        write(io, "\\section*{Аналіз централізатора (Тип пошуку: $res_type)}\n")
        
        # Вихідні дані
        l_orig_x = latexify(D_orig.polys[1], env=:raw)
        l_orig_y = latexify(D_orig.polys[2], env=:raw)
        write(io, "Початкове диференціювання:\n")
        write(io, "\\begin{equation*}\n  D = \\left( $l_orig_x \\right) \\partial_x + \\left( $l_orig_y \\right) \\partial_y\n\\end{equation*}\n\n")
        
        write(io,"Правильні розв'язки: $valid_count\n\n")
        write(io,"Неправильні розв'язки: $invalid_count\n")

        # Таблиця розподілу розмірностей
        write(io, "\\subsection*{Статистика нетривіальних розв'язків}\n")
        if !isempty(stats_table)
            sorted_k = sort(collect(keys(stats_table)))
            cols = length(sorted_k)
            write(io, "\\begin{table}[h!]\n\\centering\n\\begin{tabular}{|l|" * "c|"^cols * "}\n\\hline\n")
            write(io, "Степінь \$ k \$ & " * join(sorted_k, " & ") * " \\\\ \\hline\n")
            write(io, "Dim \$ Z(D)_k / \\mathbb{K}D \$ & " * join([stats_table[k] for k in sorted_k], " & ") * " \\\\ \\hline\n")
            write(io, "\\end{tabular}\n\\end{table}\n")
        end

        # Обрані результати
        if isempty(final_selection)
            write(io, "\\noindent\\textit{Нетривіальних розв'язків не знайдено.}\n")
        else
            status = final_selection[1].is_prop ? "пропорційний" : "НЕ ПРОПОРЦІЙНИЙ"
            write(io, "\\subsection*{Результати за критерієм \$ $res_type \$}\n")
            write(io, "Статус: \\textbf{$status}. Сумарний степінь поліномів: $(final_selection[1].sum_deg).\\\\\n")
            
            write(io, "\\begin{itemize}\n")
            for c in final_selection
                # Факторизація робить формули в LaTeX значно гарнішими
                fx = Symbolics.factorize(simplify(c.polys[1]))
                fy = Symbolics.factorize(simplify(c.polys[2]))
                write(io, "  \\item \$ D_u = ($(latexify(fx, env=:raw)))\\partial_x + ($(latexify(fy, env=:raw)))\\partial_y \$ (знайдено на \$ k = $(c.k) \$)\n")
            end
            write(io, "\\end{itemize}\n")
        end
        write(io, "\\hrule\n")
    end
    println(">>> Reports: LaTeX звіт '$filename' готовий ($res_type).")
end

end # module