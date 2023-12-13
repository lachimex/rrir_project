using LinearAlgebra
using QuadGK

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
    global a = -1
    global b = 2

    for i in 1:length(nodes)
        # Replace the following line with your function evaluation at nodes[i]
        global result = result + f((a+b)/2 + ((b-a) * nodes[i]) / 2) * weights[i]
    end

    println("Wynik: ", ((b-a)/2) * result)
    return result
end

