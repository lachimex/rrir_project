using LinearAlgebra
using QuadGK

function f(x)
    return x * x
end

function gauss_legendre_quadrature(degree)
    nodes, weights = gauss(Float64, degree)
    return nodes, weights
end

degree = 4  # Stopień kwadratury Gaussa-Legendre'a
nodes, weights = gauss_legendre_quadrature(degree)

# Wyświetlanie węzłów i wag
println("Węzły: ", nodes')
println("Wagi: ", weights')

global result = 0.0

for i in 1:length(nodes)
    # Replace the following line with your function evaluation at nodes[i]
    global result = result + f(nodes[i]) * weights[i]
end

println("Wynik: ", result)

