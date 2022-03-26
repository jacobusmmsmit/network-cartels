using Graphs
using SimpleWeightedGraphs
using StatsBase
using IterTools
using Distributions

"""
    generate_hyperbolic_graph(N::Int, α, ν, T)
Arguments:
* N: Number of nodes on the graph
* α: Related to power law P(x) = x^(-γ) by γ = 2α + 1
* ν: Dictates the total radius containing the points
* T: Temperature parameter for probability of an edge
"""
function generate_hyperbolic_graph(N::Int, α, ν, T)
    if α < 0.5
        throw(ArgumentError("α must be greater than 1/2"))
    end
    V = [(0.0, 0.0) for _ in 1:N]
    E = Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}[]
    R = 2log(N / ν)
    sources = Int[]
    dests = Int[]
    weights = Float64[]
    for i in 1:N
        θ = rand(Uniform(0, 2π))
        # r = asinh((rand()/α)*(cosh(α*R)-1))/α
        r = acosh(1 + (cosh(α * R) - 1) * rand()) / α
        V[i] = (r, θ)
    end
    for i in 1:N, j in i+1:N
        dH = disth(V[i], V[j])
        Pd = 1 / (1 + exp((1 / (2T)) * (dH - R)))
        if rand() < Pd
            push!(E, (V[i], V[j]))
            push!(sources, i)
            push!(dests, j)
            push!(weights, dH)
        end
        if rand() < Pd
            push!(E, (V[j], V[i]))
            push!(sources, j)
            push!(dests, i)
            push!(weights, dH)
        end
    end
    # Return graph and embedding
    return SimpleWeightedDiGraph(sources, dests, weights), (V, E)
end