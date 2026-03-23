using Distributed

# Додаємо воркерів (наприклад, 4 ядра)
if nprocs() == 1
    addprocs(2) 
end

@everywhere begin
    using Pkg
    Pkg.activate(".")
    using LieDerivations
    using .LieDerivations.ParallelSearch
    using Symbolics
end

# 1. Визначаємо функцію-генератор (наприклад, випадкові мономи)
@everywhere function my_generator()
    @variables x y
    vars = [x, y]
    # Випадкові степені від 0 до 5
    powers = rand(0:5, 4) 
    coeffs = rand(-10:10, 2)
    
    P1 = coeffs[1] * x^powers[1] * y^powers[2]
    P2 = coeffs[2] * x^powers[3] * y^powers[4]
    
    return Derivation([P1, P2], vars)
end

println("Початок паралельного аналізу на $(nworkers()) воркерах...")

# 2. Запускаємо масовий тест
# Знайдемо централізатори для 50 випадкових диференціювань до степеня k=5
t_start = time()

# Запуск масового тесту
all_results = ParallelSearch.run_parallel_tests(my_generator, 10, 10)

total_exec_time = time() - t_start

# 3. Аналіз результатів
println("\nЗнайдено нетривіальних рішень:")
for (D, found) in all_results
    interesting = filter(f -> f.is_interesting, found)
    if !isempty(interesting)
        println("Для D = $(D.polys) знайдено $(length(interesting)) цікавих розв'язків.")
    end
end



# Генеруємо фінальний звіт
println("Формування звіту...")
Reports.generate_mass_report(all_results, "tex-code/mass_report.tex", total_exec_time)

println("Готово. Звіт збережено в 'tex-code/mass_report.tex'")