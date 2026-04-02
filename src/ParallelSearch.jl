module ParallelSearch

using ..LieDerivations
using ThreadsX # Потрібно додати через Pkg.add("ThreadsX")

export run_threaded_tests

"""
Функція для багатопотокового пошуку централізаторів (економія RAM).
"""
function run_threaded_tests(generator_func::Function, total_tests::Int, max_k::Int)
    println(">>> Запуск потокового пошуку: $total_tests тестів на $(Threads.nthreads()) потоках.")
    
    tasks = 1:total_tests

    # ThreadsX.map працює аналогічно до pmap, але в межах одного процесу
    results_list = ThreadsX.map(tasks) do i
        D_rand = generator_func()
        # Викликаємо ваш стандартний solve_centralizer
        search_res = solve_centralizer(D_rand, max_k)
        
        if i % 10 == 0
            println("Thread $(Threads.threadid()): Оброблено $i тестів...")
        end

        # Повертаємо пару для формування словника
        return D_rand => search_res.all_res
    end

    # Перетворюємо масив пар у словник
    return Dict(results_list)
end

end # module