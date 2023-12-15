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

function B(dw_dx, dv_dx)
    result, error = quadgk(x -> dw_dx(x) * dv_dx(x), 0, 3)
    return -1 * result
end

function e_i_v1(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 2)
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
    x_1 = dx * div(i-1, 2)
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

function de_dx(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 2)
    x_2 = x_1 + dx
    a = 1/(x_2 - x_1)
    return function (x)
        if x_1 <= x <= x_2
            return a
        else
            return 0
        end
    end
end
    

n = 7
A = zeros(Float64, 2*n, 2*n)
C = zeros(Float64, 2*n, 1)
D = Array{Function}(undef, 1, 2*n)
E = Array{Function}(undef, 1, 2*n)

for i in 1:2*n
    if i % 2 == 0
        e_i = e_i_v1(n, i)
    else
        e_i = e_i_v2(n, i)
    end
    E[i] = e_i
    D[i] = de_dx(n, i)
end

for i in 1:2*n
    for j in i:2*n
        if (j - i) < 2
            A[i, j] = B(D[i], D[j])
            A[j, i] = A[i, j]
        end
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
    for i in 1:2*n
        out += W[i] * E[i](x)
    end
    return out
end

function u(x)
    out = 5 - x/3
    out += w(x)
    return out
end

x = range(0.75, 2.25, length=1000)
y = u.(x)
gui(plot(x, y))
sleep(1000000)