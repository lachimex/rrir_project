using LinearAlgebra
using QuadGK
using ForwardDiff
using Plots

const G = 6.67430e-11

function ro(x)
    if (1 < x && x <= 2)
        return 1
    else
        return 0
    end
end

function L(v)
    result, error = quadgk(x -> v(x) * ro(x), 0, 3)
    return 4 * pi * G * result;
end

function B(w, v)
    result, error = quadgk(x -> ForwardDiff.derivative(w, x) * ForwardDiff.derivative(v, x), 0, 3)
    return -1 * result;
end

function e_i_v1(n, i)
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

function e_i_v2(n, i)
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
    if i % 3 == 1
        e_i = e_i_v1(n, i)
    elseif i % 3 == 2
        e_i = e_i_v2(n, i)
    else
        e_i = e_quadratic(n, i)
    end
    E[i] = e_i
end

for i in 1:3*n
    for j in i:3*n
        A[i, j] = B(E[i], E[j])
        A[j, i] = A[i, j]
    end
    C[i] = L(E[i])
end

println()
for i in 1:size(C, 1)
    println(C[i, :])
end

println()
for i in 1:size(A, 1)
    println(A[i, :])
end
println()

W = C \ A
println(W)

function w(x)
    out = 0.0
    for i in 1:3*n
        out += W[i] * E[i](x)
    end
    return out
end

function u(x)
    out = 5 - x/3
    out += w(x)
    return out
end

println(w(0))
println(w(3))
x = range(0, 3, length=1000)
y = u.(x)
gui(plot(x, y))
sleep(1000000)