using LinearAlgebra
using QuadGK
using ForwardDiff
using Plots

const G = 6.67430e-11

function f(x)
    return x * x
end

function ro(x)
    if (1 < x && x <= 2)
        return 1
    else
        return 0
    end
end

function integral_on_a_b(fun, a, b)
    degree = 10
    nodes, weights = gauss(Float64, degree)

    global result = 0.0
    for i in 1:length(nodes)
        x = (a+b)/2 + ((b-a) * nodes[i]) / 2
        global result = result + ro(x) * fun(x) * weights[i]
    end

    return ((b-a)/2) * result
end

function integral_on_a_b_v2(fun1, fun2, a, b)
    degree = 10
    nodes, weights = gauss(Float64, degree)

    global result = 0.0
    for i in 1:length(nodes)
        x = (a+b)/2 + ((b-a) * nodes[i]) / 2
        global result = result + ForwardDiff.derivative(fun1, x) * ForwardDiff.derivative(fun2, x) * weights[i]
    end

    return ((b-a)/2) * result
end

function L(v)
    degree = 5
    result, error = quadgk(x -> v(x) * ro(x), 0, 3)
    return 4 * pi * G * result;
end

function B(w, v)
    degree = 5
    result, error = quadgk(x -> ForwardDiff.derivative(w, x) * ForwardDiff.derivative(v, x), 0, 3)
    return -1 * result;
end

function f1(x)
    return x
end

function e_i_odd(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 3)
    x_2 = x_1 + dx
    a = 1/(x_1 - x_2)
    b = 1 - x_1/(x_1 - x_2)
    # println("e$i $a x + $b dla przedzialu ($x_1, $x_2)")
    return function (x)
        if x_1 <= x <= x_2
            return a * x + b
        else
            return 0
        end
    end
end

function e_i_even(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 3)
    x_2 = x_1 + dx
    a = 1/(x_2 - x_1)
    b = 1 - x_2/(x_2 - x_1)
    # println("e$i $a x + $b dla przedzialu ($x_1, $x_2)")
    return function (x)
        if x_1 <= x <= x_2
            return a * x + b
        else
            return 0
        end
    end
end

function e_quadratic(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 3)
    x_2 = x_1 + dx
    return function (x)
        if x_1 <= x <= x_2
            return -1 * (x-x_1) * (x-x_2)
        else
            return 0
        end
    end
end

n = 3
A = zeros(Float64, 3*n, 3*n)
C = zeros(Float64, 3*n, 1)
E = Array{Function}(undef, 1, 3*n)

for i in 1:3*n
    if i % 3 == 0
        e_i = e_i_even(n, i)
    elseif i % 3 == 1
        e_i = e_i_odd(n, i)
    else
        e_i = e_quadratic(n, i)
    end
    E[i] = e_i
end

for i in 1:3*n
    for j in i:3*n
        A[i, j] = B(E[i], E[j])
    end
    C[i] = L(E[i])
end

println()
for i in 1:size(C, 1)
    println(C[i, :])
end

println()
A = 0.5 * (A + A')
for i in 1:size(A, 1)
    println(A[i, :])
end
println()

W = C \ A
println(W)

function u(x)
    out = 5 - x/3
    for i in 1:3*n
        out += W[i] * E[i](x)
    end
    return out
end

x = range(1, 2, length=1000)
y = u.(x)
gui(plot(x, y))
sleep(1000000)