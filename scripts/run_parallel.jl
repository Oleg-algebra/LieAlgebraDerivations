using Pkg
Pkg.activate(".")

using LieDerivations
using .LieDerivations.ParallelSearch
using Symbolics

# 1. Визначаємо функцію-генератор (тепер просто звичайна функція)
function my_generator()
    @variables x y
    vars = [x, y]
    powers = rand(0:5, 4) 
    coeffs = rand(-10:10, 2)
    
    P1 = coeffs[1] * x^powers[1] * y^powers[2]
    P2 = coeffs[2] * x^powers[3] * y^powers[4]
    
    return Derivation([P1, P2], vars)
end

# Налаштування параметрів
total_number_of_tests = 100
max_degree_K = 5

println("Початок потокового аналізу на $(Threads.nthreads()) потоках...")

# 2. Запуск масового тесту з заміром часу
t_start = time()

# Використовуємо нову функцію run_threaded_tests
all_results = ParallelSearch.run_threaded_tests(my_generator, total_number_of_tests, max_degree_K)

total_exec_time = time() - t_start

# 3. Швидкий аналіз результатів у консолі
println("\nЗнайдено нетривіальних рішень:")
for (D, found) in all_results
    interesting = filter(f -> f.is_interesting, found)
    if !isempty(interesting)
        println("Для D = $(D.polys) знайдено $(length(interesting)) цікавих розв'язків.")
    end
end

# 4. Формування LaTeX звіту
println("\nФормування звіту...")
Reports.generate_mass_report(all_results, "tex-code/mass_report.tex", total_exec_time)

println("Готово. Звіт збережено в 'tex-code/mass_report.tex'")