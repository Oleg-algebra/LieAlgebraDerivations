using LieDerivations
using Symbolics
using Latexify
using LinearAlgebra
using Revise
using LaTeXStrings # пакет для роботи з сирим LaTeX кодом
using Dates       

# 1. Ініціалізація змінних
@variables x y

# 2. Визначення вашого PhD диференціювання (приклад Jacobian)
# Нехай D = (x^2)*∂/∂x + (2xy)*∂/∂y
poly_x = y
poly_y = x
D = Derivation([poly_x, poly_y], [x, y])

max_deg = 15
println("Шукаємо централізатор для D до степеня $max_deg...")

# # Перший раз (з компіляцією)
# @time ns, coeffs, u_polys = solve_centralizer(D, max_deg)

# # Другий раз (чиста швидкість)
@time ns, coeffs, u_polys = solve_centralizer(D, 1)

t_exec = @elapsed begin
    global ns, coeffs, u_polys = solve_centralizer(D, max_deg)
end
println("Час виконання: $(round(t_exec, digits=3)) сек.")

# 4. Форматування результату в LaTeX
if isempty(ns)
    println("Централізатор складається тільки з тривіальних розв'язків.")
else
    println("\nЗнайдено базисні елементи централізатора:")
    
    for col in 1:size(ns, 2)
        # 1. Підставляємо знайдені коефіцієнти
        sol_dict = Dict(coeffs .=> ns[:, col])
        final_u_polys = [substitute(p, sol_dict) for p in u_polys]
        
        # 2. Перетворюємо поліноми в LaTeX-рядки окремо
        # Використовуємо env=:raw, щоб отримати чистий LaTeX без оточення $...$
        l_poly_x = latexify(final_u_polys[1], env=:raw)
        l_poly_y = latexify(final_u_polys[2], env=:raw)
        
        # 3. Формуємо фінальний рядок, використовуючи макрос L"..." (LatexString)
        # Це дозволяє вставляти змінні прямо в LaTeX код
        final_latex = L"D_{u_{%$col}} = (%$l_poly_x) \frac{\partial}{\partial x} + (%$l_poly_y) \frac{\partial}{\partial y}"
        
        println("--- Базисний елемент $col ---")
        display(final_latex) # VS Code відрендерить це як гарну формулу
        println("\nСирий код для LaTeX:")
        println(final_latex)
    end
end

filename = "tex-code/results.tex"

open(filename, "w") do io
    # 1. Службова інформація
    write(io, "% !TeX root = main.tex\n")
    # Додамо час виконання у LaTeX звіт
    write(io, "\\noindent\\textit{Час обчислення: $(round(t_exec, digits=3)) сек.}\\\\\n")
    write(io, "% Розрахунок виконано: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n\n")
    
    # 2. Опис початкового диференціювання D
    # Перетворюємо поліноми D у LaTeX
    l_orig_x = latexify(D.polys[1], env=:raw)
    l_orig_y = latexify(D.polys[2], env=:raw)
    
    write(io, "\\section*{Централізатор для диференціювання}\n")
    write(io, "\\begin{equation*}\n")
    write(io, "  D = \\left( $l_orig_x \\right) \\frac{\\partial}{\\partial x} + \\left( $l_orig_y \\right) \\frac{\\partial}{\\partial y}\n")
    write(io, "\\end{equation*}\n")
    
    write(io, "Шукаємо централізатор у просторі поліномів до степеня \$k = $max_deg\$.\n")
    write(io, "\\vspace{1em}\n\n")

    # 3. Вивід знайдених базисних елементів
    if isempty(ns)
        write(io, "Централізатор у заданому просторі містить лише тривіальні розв'язки.\n")
    else
        write(io, "\\textbf{Базис централізатора:}\n")
        write(io, "\\begin{itemize}\n")
        
        for col in 1:size(ns, 2)
            # Підстановка та формування кінцевого полінома
            sol_dict = Dict(coeffs .=> ns[:, col])
            final_u_polys = [substitute(p, sol_dict) for p in u_polys]

            factored_x = Symbolics.factorize(final_u_polys[1])
            factored_y = Symbolics.factorize(final_u_polys[2])
            
            # Очищення від дуже малих значень або спрощення
            l_poly_x = latexify(factored_x, env=:raw)
            l_poly_y = latexify(factored_y, env=:raw)
            
            latex_entry = "  \\item \$ D_{u_{$col}} = \\left( $l_poly_x \\right) \\frac{\\partial}{\\partial x} + \\left( $l_poly_y \\right) \\frac{\\partial}{\\partial y} \$\n"
            write(io, latex_entry)
        end
        
        write(io, "\\end{itemize}\n")
    end
    
    write(io, "\\hrule\\vspace{2em}\n\n")
end

println("Звіт для PhD успішно згенеровано у $filename")