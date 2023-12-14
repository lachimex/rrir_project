using LinearAlgebra
using QuadGK
using ForwardDiff

const G = 6.67430e-11
const h = 1e-8

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
    degree = 50
    nodes, weights = gauss(Float64, degree)

    global result = 0.0
    for i in 1:length(nodes)
        x = (a+b)/2 + ((b-a) * nodes[i]) / 2
        global result = result + ro(x) * fun(x) * weights[i]
    end

    return ((b-a)/2) * result
end

function integral_on_a_b_v2(fun1, fun2, a, b)
    degree = 50
    nodes, weights = gauss(Float64, degree)

    global result = 0.0
    for i in 1:length(nodes)
        x = (a+b)/2 + ((b-a) * nodes[i]) / 2
        global result = result + ForwardDiff.derivative(fun1, x) * ForwardDiff.derivative(fun2, x) * weights[i]
    end

    return ((b-a)/2) * result
end

function L(v)
    return integral_on_a_b(v, 0, 3);
end

function B(w, v)
    return -1 * integral_on_a_b_v2(w, v, 0, 3);
end

n = 2

function e_i_odd(n, i)
    dx = 3 / n
    x_1 = dx * div(i, 2)
    x_2 = x_1 + dx
    a = 1/(x_1 - x_2)
    b = 1 - x_1/(x_1 - x_2)
    return function (x)
        a * x + b
    end
end

function e_i_even(n, i)
    dx = 3 / n
    x_1 = dx * div(i, 2)
    x_2 = x_1 + dx
    a = 1/(x_2 - x_1)
    b = 1 - x_2/(x_2 - x_1)
    return function (x)
        a * x + b
    end
end

A = zeros(Float64, 2*n, 2*n)
C = zeros(Float64, 2*n, 1)

for i in 1:2*n
    if i % 2 == 0
        e_i = e_i_even(n, i)
    else
        e_i = e_i_odd(n, i)
    end
    for j in 1:2*n
        if j % 2 == 0
            e_j = e_i_even(n, j)
        else
            e_j = e_i_odd(n, j)
        end
        A[i, j] = B(e_i, e_j)
    end
    C[i] = L(e_i)
end

println(C \ A)