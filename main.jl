using LinearAlgebra
using QuadGK
using ForwardDiff

const G = 6.67430e-11
const h = 1e-8

function f(x)
    return x * x
end

function integral_on_a_b(fun, a, b)
    degree = 10
    nodes, weights = gauss(Float64, degree)

    global result = 0.0
    for i in 1:length(nodes)
        global result = result + fun((a+b)/2 + ((b-a) * nodes[i]) / 2) * weights[i]
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
    return 4 * π * G * integral_on_a_b(v, 1, 2);
end

function B(w, v)
    return -1 * integral_on_a_b_v2(w, v, 0, 3);
end

function f1(x)
    return x * x
end

function f2(x)
    return x * x
end

println(B(f1, f2))

