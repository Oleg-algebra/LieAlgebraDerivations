# 1. Налаштування середовища
using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) # Активація проекту (коренева папка)

using Revise           # Автоматичне оновлення модулів при зміні коду
using Symbolics
using LieDerivations
using .LieDerivations.Reports # Доступ до модуля Reports через головний модуль

# 2. Визначення задачі (Приклад)
println("--- Ініціалізація обчислень ---")
@variables x y
vars = [x, y]

# Вставте тут ваше диференціювання D
# Приклад: D = x^2 * ∂/∂x + y^2 * ∂/∂y
D_orig = Derivation([y^2,x^2 ], vars) 

max_deg = 10  # Максимальний степінь для пошуку
filename = "tex-code/centralizer_report.tex" # Шлях до файлу звіту

# 3. Запуск основного солвера
println("Пошук централізатора до степеня k = $max_deg...")
t_start = time()

# solve_centralizer тепер повертає об'єкт з полем .all_found (список усіх розв'язків)
res_data = solve_centralizer(D_orig, max_deg)

t_exec = time() - t_start
println("Обчислення завершено за $(round(t_exec, digits=2)) сек.")

# 4. Генерація звіту через модуль Reports
# Можна обрати res_type = "min" або "max"
println("Генерація LaTeX звіту ($filename)...")

if !isempty(res_data.all_res)
    Reports.generate_latex_report(
        res_data.all_res, 
        D_orig, 
        vars, 
        filename; 
        res_type = "all"
    )
    println(">>> Звіт успішно збережено.")
else
    println("!!! Нетривіальних розв'язків не знайдено, звіт не згенеровано.")
end

println("--- Роботу завершено ---")