using LinearAlgebra
using QuadGK
using Math
using ForwardDiff

const G = 6.67430e-11

function f(x)
    return x * x
end

function gauss_legendre_quadrature(degree)
    nodes, weights = gauss(Float64, degree)
    return nodes, weights
end

function integral_on_a_b(fun, a, b)
    degree = 10  # Stopień kwadratury Gaussa-Legendre'a
    nodes, weights = gauss_legendre_quadrature(degree)

    # Wyświetlanie węzłów i wag
    println("Węzły: ", nodes')
    println("Wagi: ", weights')

    global result = 0.0
    for i in 1:length(nodes)
        # Replace the following line with your function evaluation at nodes[i]
        global result = result + fun((a+b)/2 + ((b-a) * nodes[i]) / 2) * weights[i]
    end

    return ((b-a)/2) * result
end

function L(v)
    return 4 * pi * G * integral_on_a_b(v, 1, 2);
end

function B(w, v)
    return -1 * integral_on_a_b(ForwardDiff.derivative(w) * ForwardDiff.derivative(v), 0, 3);
end

