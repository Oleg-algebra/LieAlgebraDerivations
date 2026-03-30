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
        if res_type == "all"
            # Повертаємо всі знайдені кандидати без врахування степеня
            final_selection = candidates
        else
            # Стара логіка для відбору за екстремальним степенем
            all_degs = [c.sum_deg for c in candidates]
            target = (res_type == "min") ? minimum(all_degs) : maximum(all_degs)
            final_selection = filter(c -> c.sum_deg == target, candidates)
        end
    end

    # 3. Запис у LaTeX
    open(filename, "w") do io
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



"""
Генерує загальний звіт для масового паралельного тестування.
Аргументи:
- results: Dict{Derivation, Vector} (результат ParallelSearch)
- filename: назва .tex файлу
- total_time: загальний час виконання (с)
"""
function generate_mass_report(results, filename::String, total_time::Float64)
    open(filename, "w") do io
        # --- ПРЕАМБУЛА ТА ШАПКА ---
        write(io, "\\section*{Масовий аналіз централізаторів Lie-диференціювань}\n")
        
        # Розрахунок загальної статистики
        total_tests = length(results)
        unprop_count = 0
        inv_count = 0
        for (_, found) in results
            unprop_count += count(s -> s.is_interesting, found)
            inv_count += count(s -> get(s, :is_invariant, false), found)
        end

        write(io, "\\begin{tcolorbox}[colback=blue!5,colframe=blue!75,title=Метрики обчислень]\n")
        write(io, "\\begin{tabular}{ll}\n")
        write(io, "\\textbf{Кількість тестів:} & $total_tests \\\\\n")
        write(io, "\\textbf{Знайдено комутаторів (\$D_u \\notin \\mathbb{K}[x,y]D\$):} & \\textbf{$unprop_count} \\\\\n")
        write(io, "\\textbf{Знайдено інваріантів (\$f \\cdot D, D(f)=0\$):} & $inv_count \\\\\n")
        write(io, "\\hline\n")
        write(io, "\\textbf{Загальний час:} & $(round(total_time, digits=2)) с \\\\\n")
        write(io, "\\textbf{Сер. час на тест:} & $(round(total_time/total_tests, digits=4)) с \\\\\n")
        write(io, "\\end{tabular}\n")
        write(io, "\\end{tcolorbox}\n\n")

        # --- ДЕТАЛЬНИЙ ВИВІД ---
        for (D, found_list) in results
            l_orig = "\$ D = (" * latexify(D.polys[1], env=:raw) * ")\\partial_x + (" * latexify(D.polys[2], env=:raw) * ")\\partial_y \$"
            
            write(io, "\\subsection*{Аналіз поля: $l_orig}\n")
            write(io, "\\begin{itemize}[leftmargin=1cm]\n")

            # 1. Секція Комутаторів (is_interesting)
            interesting = filter(s -> s.is_interesting, found_list)
            if !isempty(interesting)
                write(io, "  \\item \\textbf{Централізатор (непропорційні):}\n")
                write(io, "    \\begin{itemize}\n")
                for s in interesting
                    write(io, "      \\item \$ D_u = ($(latexify(s.polys[1], env=:raw)))\\partial_x + ($(latexify(s.polys[2], env=:raw)))\\partial_y \$ (k=$(s.degree))\n")
                end
                write(io, "    \\end{itemize}\n")
            end

            # 2. Секція Інваріантних кратних (is_invariant)
            invariants = filter(s -> get(s, :is_invariant, false), found_list)
            if !isempty(invariants)
                write(io, "  \\item \\textbf{Інваріантні кратні (\$f \\cdot D, D(f)=0\$):}\n")
                write(io, "    \\begin{itemize}\n")
                for s in invariants
                    write(io, "      \\item \$ D_u = ($(latexify(s.polys[1], env=:raw)))\\partial_x + ($(latexify(s.polys[2], env=:raw)))\\partial_y \$ (k=$(s.degree))\n")
                end
                write(io, "    \\end{itemize}\n")
            end

            # 3. Базовий випадок (пропорційне k=deg)
            base = filter(s -> !s.is_interesting && !get(s, :is_invariant, false), found_list)
            if !isempty(base)
                write(io, "  \\item \\textit{Базове поле \$D\$ підтверджено на k=$(base[1].degree).}\n")
            end

            write(io, "\\end{itemize}\n")
            write(io, "\\hrulefill\n")
        end
    end
end



end # module