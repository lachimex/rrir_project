using LinearAlgebra
using QuadGK
using Plots

const G = 6.67430e-11

function ro(x)
    if (1 < x && x <= 2)
        return 1
    else
        return 0
    end
end

function L(v, x_1, x_2)
    result, error = quadgk(x -> v(x) * ro(x), x_1, x_2)
    return 4 * pi * G * result;
end

function B(dw_dx, dv_dx, x_1, x_2)
    result, error = quadgk(x -> dw_dx(x) * dv_dx(x), x_1, x_2)
    return -1 * result
end

function e_i_v1(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 2)
    x_2 = x_1 + dx
    function f(x)
        if x_1 <= x <= x_2
            return x / n - i + 1
        else
            return 0
        end
    end
    return f, x_1, x_2
end

function e_i_v2(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 2)
    x_2 = x_1 + dx
    function f(x)
        if x_1 <= x <= x_2
            return -x / n + + 1
        else
            return 0
        end
    end
    return f, x_1, x_2
end

function de_dx(n, i)
    dx = 3 / n
    x_1 = dx * div(i-1, 2)
    x_2 = x_1 + dx
    if i % 2 == 1
        a = 1/n
    else
        a = -1/n
    end
    function f(x)
        if x_1 <= x <= x_2
            return a
        else
            return 0
        end
    end
    return f, x_1, x_2
end
    

n = 15
A = zeros(Float64, 2*n, 2*n)
C = zeros(Float64, 2*n, 1)
D = Array{Function}(undef, 2*n, 1) #array of derivatives
D2 = zeros(Float64, 2*n, 2) #array of boundaries
E = Array{Function}(undef, 2*n, 1) #array of functions

for i in 1:2*n
    if i % 2 == 0
        e_i, x_1, x_2 = e_i_v1(n, i)
    else
        e_i, x_1, x_2 = e_i_v2(n, i)
    end
    E[i] = e_i
    d, x_1, x_2 = de_dx(n, i)
    D[i] = d
    D2[i, 1] = x_1
    D2[i, 2] = x_2
end

for i in 1:2*n
    for j in i:2*n
        if (j - i) < 2
            A[i, j] = B(D[i], D[j], D2[i, 1], D2[i, 2])
            A[j, i] = A[i, j]
        end
    end
    C[i] = L(E[i], D2[i, 1], D2[i, 2])
end


println()
for i in 1:size(A, 1)
    println(A[i, :])
end
println()
for i in 1:size(C, 1)
    println(C[i, :])
end
println()

W = pinv(A) * C
for i in 1:size(W, 1)
    println(W[i, :])
end

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

x = range(0, 3, length=1000)
y = u.(x)
gui(plot(x, y))
sleep(100)
