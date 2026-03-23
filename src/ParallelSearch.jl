module ParallelSearch

using Distributed
using ..LieDerivations # Звернення до батьківського модуля
using Symbolics

export run_parallel_tests

"""
Функція для паралельного пошуку централізаторів.
Аргументи:
- generator_func: функція без аргументів, що повертає об'єкт Derivation
- total_tests: загальна кількість диференціювань, які треба дослідити
- max_k: максимальний степінь пошуку для кожного теста
"""
function run_parallel_tests(generator_func::Function, total_tests::Int, max_k::Int)
    println(">>> Запуск паралельного пошуку: $total_tests тестів на $(nworkers()) воркерах.")
    
    # Створюємо масив завдань (просто індекси)
    tasks = 1:total_tests

    # pmap автоматично розподіляє завдання між воркерами
    # Кожен воркер виконує анонімну функцію
    results_list = pmap(tasks) do i
        # 1. Генеруємо випадкове диференціювання
        D_rand = generator_func()
        
        # 2. Шукаємо централізатор (solve_centralizer має повертати NamedTuple з all_found)
        # Припускаємо, що solve_centralizer повертає (all_found=..., is_valid=...)
        search_res = solve_centralizer(D_rand, max_k)
        
        # 3. Повертаємо пару: Диференціювання => Його знайдені комутатори
        return D_rand => search_res.all_res
    end

    # Перетворюємо список пар у словник для зручного доступу
    return Dict(results_list)
end

end # module