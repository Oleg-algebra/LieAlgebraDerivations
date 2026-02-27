using LieDerivations
using Symbolics
using BenchmarkTools

# 1. Визначаємо змінні
@variables x y

# 2. Створюємо приклад диференціювання (наприклад, якобієве)
# D = y*∂/∂x + 0*∂/∂y
D = Derivation([y, 0], [x, y])

# 3. Визначаємо потенційний комутатор Du
Du = Derivation([x, y + x^2], [x, y])

lie_bracket(D,Du)