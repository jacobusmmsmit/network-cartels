# include("../src/NetworkCartelDetection.jl")
# using .NetworkCartelDetection

# Fundamental data structures
using Graphs, SimpleWeightedGraphs
# Reading files
using GraphIO, DelimitedFiles
# Others
using StatsBase
using Plots
using IterTools
using Distributions
using Measures

using Random
using SparseArrays

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
    # return SimpleWeightedDiGraph(sources, dests, weights), (V, E)

    return SimpleDiGraph([Edge(i) for i in zip(sources, dests)]), (V, E)
end

disth((r, θ), (r′, θ′)) = acosh(cosh(r) * cosh(r′) - sinh(r) * sinh(r′) * cos(π - abs(π - abs(θ - θ′))))

function get_dm(vertices, distance_function=diste)
    [p1 == p2 ? 0.0 : distance_function(p1, p2) for (p1, p2) in Iterators.product(vertices, vertices)]
end

function _p2p_cost(i, strategy, sp_dm, geo_dm, α)
    cost = α * length(strategy[strategy.>0])
    for j in Iterators.flatten((1:(i-1), (i+1):length(strategy)))
        cost += sp_dm[i, j] / geo_dm[i, j]
    end
    return cost
end

function evaluate_cost(cost_function, adj_matrix)
    map(i_s -> cost_function(i_s[1], i_s[2]), enumerate(eachrow(adj_matrix)))
end

function social_cost(cost_function, adj_matrix)
    mapreduce(i_s -> cost_function(i_s[1], i_s[2]), +, enumerate(eachrow(adj_matrix)))
end

function calculate_p2p_cost(graph, embedded_distance_matrix, edge_cost, ds=floyd_warshall_shortest_paths(graph, embedded_distance_matrix), shortestpath_distance_matrix=ds.dists)
    cfn = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
    evaluate_cost(cfn, adjacency_matrix(graph))
end

function calculate_p2p_social_cost(graph, embedded_distance_matrix, edge_cost, ds=floyd_warshall_shortest_paths(graph, embedded_distance_matrix), shortestpath_distance_matrix=ds.dists)
    cfn = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
    social_cost(cfn, adjacency_matrix(graph))
end

function best_response(player, _adj_matrix, _shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
    adj_matrix = copy(_adj_matrix)
    shortestpath_distance_matrix = copy(_shortestpath_distance_matrix)
    cost_function = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
    min_cost = evaluate_cost(cost_function, adj_matrix)[player]
    BR = adj_matrix[player, :]
    i = 1
    for _strategy in subsets(collect(Iterators.flatten((1:(player-1), (player+1):size(adj_matrix, 1)))))
        strategy = Vector(sparsevec(_strategy, 1, size(adj_matrix, 1)))
        adj_matrix[player, :] = strategy
        dropzeros!(adj_matrix)
        g = SimpleDiGraph(adj_matrix)
        ds = floyd_warshall_shortest_paths(g, embedded_distance_matrix)
        shortestpath_distance_matrix = ds.dists
        cost_function = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
        cur_cost = evaluate_cost(cost_function, adj_matrix)[player]
        if cur_cost < min_cost
            min_cost = cur_cost
            BR = strategy
        end
        i += 1
    end
    return BR, min_cost
end


function best_response(player, cartel, _adj_matrix, _shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
    adj_matrix = copy(_adj_matrix)
    shortestpath_distance_matrix = copy(_shortestpath_distance_matrix)
    cost_function = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
    min_cost = mapreduce(player -> evaluate_cost(cost_function, adj_matrix)[player], +, cartel)
    BR = adj_matrix[player, :]
    i = 1
    for _strategy in collect(subsets(collect(Iterators.flatten((1:(player-1), (player+1):size(adj_matrix, 1))))))
        strategy = Vector(sparsevec(_strategy, 1, size(adj_matrix, 1)))
        adj_matrix[player, :] = strategy
        dropzeros!(adj_matrix)
        g = SimpleDiGraph(adj_matrix)
        ds = floyd_warshall_shortest_paths(g, embedded_distance_matrix)
        shortestpath_distance_matrix = ds.dists
        cost_function = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
        cur_cost = mapreduce(player -> evaluate_cost(cost_function, adj_matrix)[player], +, cartel)
        if cur_cost < min_cost
            min_cost = cur_cost
            BR = strategy
        end
        i += 1
    end
    return BR, evaluate_cost(cost_function, adj_matrix)[player]
end

function cartel_best_response(cartel, _adj_matrix, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost; n=2)
    new_adj = copy(_adj_matrix)
    for player in cartel
        cost_function = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
        if isinf(evaluate_cost(cost_function, new_adj)[player])
            new_adj[player, :] = best_response(player, cartel, new_adj, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)[1]
        end
        if isinf(evaluate_cost(cost_function, new_adj)[player])
            new_adj[player, :] .= 1
        end
    end
    for _ in 1:n
        for player in cartel
            new_adj[player, :] = best_response(player, cartel, new_adj, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)[1]
        end
    end
    return new_adj
end

function best_response_dynamics(adj_matrix, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost; n=1, cartels=Vector{Vector{Int}}[])
    new_adj = copy(adj_matrix)
    if length(cartels) == 0
        for _ in 1:n
            for player in 1:size(adj_matrix, 1)
                new_adj[player, :] = best_response(player, new_adj, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)[1]
                dropzeros!(new_adj)
            end
        end
    else
        cartel_nodes = reduce(vcat, cartels)
        normal_nodes = setdiff(1:size(adj_matrix, 1), cartel_nodes)
        for _ in 1:n
            for player in normal_nodes
                new_adj[player, :] = best_response(player, new_adj, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)[1]
                dropzeros!(new_adj)
            end
            for cartel in cartels
                new_adj = cartel_best_response(cartel, new_adj, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
            end
        end
    end
    return new_adj
end

"""
    nz(A, B)
Percentage of differences in two inputs `A` and `B`
"""
nz(A, B) = sum(A .!= B) / length(A)

# Please accept my apologies, thou who dost gaze upon this code
function analyse_links(A, B, cartel, by=:hamming)
    nocartel = setdiff(1:size(A, 1), cartel)
    valid = true
    if nocartel != collect(5:10)
        valid = false
    end
    if by == :hamming
        cc = nz(A[cartel, cartel], B[cartel, cartel])
        cnc = nz(A[cartel, nocartel], B[cartel, nocartel])
        ncc = nz(A[nocartel, cartel], B[nocartel, cartel])
        ncnc = nz(A[nocartel, nocartel], B[nocartel, nocartel])
    elseif by == :nconnections
        cc = sum(A[cartel, cartel]) / sum(B[cartel, cartel])
        cnc = sum(A[cartel, nocartel]) / sum(B[cartel, nocartel])
        ncc = sum(A[nocartel, cartel]) / sum(B[nocartel, cartel])
        ncnc = sum(A[nocartel, nocartel]) / sum(B[nocartel, nocartel])
    elseif by == :An
        cc = sum(A[cartel, cartel])
        cnc = sum(A[cartel, nocartel])
        ncc = sum(A[nocartel, cartel])
        ncnc = sum(A[nocartel, nocartel])
    elseif by == :Bn
        cc = sum(B[cartel, cartel])
        cnc = sum(B[cartel, nocartel])
        ncc = sum(B[nocartel, cartel])
        ncnc = sum(B[nocartel, nocartel])
    end
    res = [cc cnc
        ncc ncnc]
    # map!(x -> isnan(x) ? zero(x) : x, res, res)
    return res, valid
end

function main()
    Random.seed!(2)
    # toy = loadgraph("data/toy_network.txt", "graph_key", EdgeListFormat())
    # toy_embedding = readdlm("data/toy_network_polar.csv", ',', skipstart=1, Float64)[:, 2:3] |>
    #                 x -> [(r, theta) for (r, theta) in eachrow(x)]

    function apply_method()
        N = 10
        toy, (toy_embedding, _) = generate_hyperbolic_graph(N, 0.7, 1.0, 0.15)
        toy_dm = get_dm(toy_embedding, disth)
        graph = toy
        embedded_distance_matrix = toy_dm
        edge_cost = 1.0
        ds = floyd_warshall_shortest_paths(graph, embedded_distance_matrix)
        shortestpath_distance_matrix = ds.dists
        # cfn = (player, strategy) -> _p2p_cost(player, strategy, shortestpath_distance_matrix, embedded_distance_matrix, edge_cost)
        new_adj = best_response_dynamics(adjacency_matrix(graph), shortestpath_distance_matrix, embedded_distance_matrix, edge_cost, n=2)
        new_adj_cartels = best_response_dynamics(adjacency_matrix(graph), shortestpath_distance_matrix, embedded_distance_matrix, edge_cost, n=2, cartels=[[1, 2, 3, 4]])

        # new_adj_opt = best_response_dynamics(adjacency_matrix(graph), shortestpath_distance_matrix, embedded_distance_matrix, edge_cost, n=6, cartels=[collect(1:9)])
        # adjacency_matrix(graph)
        # display(calculate_p2p_social_cost(SimpleDiGraph(new_adj_cartels), embedded_distance_matrix, edge_cost))
        # display(calculate_p2p_social_cost(SimpleDiGraph(new_adj_opt), embedded_distance_matrix, edge_cost))
        # display(new_adj)
        # p1 = plot_embedded(toy, toy_embedding)
        # p2 = plot_embedded(SimpleDiGraph(new_adj), toy_embedding)
        # p3 = plot_embedded(SimpleDiGraph(new_adj_cartels), toy_embedding)
        # p4 = plot_embedded(SimpleDiGraph(new_adj_opt), toy_embedding)
        # plot(p1, p2, p3, p4, size=(600, 600))
        # # savegraph("out/toy_BR_network.txt", SimpleDiGraph(new_adj))
        # old_adj = adjacency_matrix(graph)
        suspiciousness = zeros(Float64, length(toy.fadjlist))
        for (i, (old_strat, new_strat)) in enumerate(zip(eachrow(new_adj), eachrow(new_adj_cartels)))
            suspiciousness[i] = sum(old_strat .!= new_strat)
        end
        suspiciousness .= suspiciousness ./ sum(suspiciousness)
        ham, vdham = analyse_links(new_adj, new_adj_cartels, 1:4, :hamming)
        nc, vdnc = analyse_links(new_adj, new_adj_cartels, 1:4, :nconnections)
        An, vdB = analyse_links(new_adj, new_adj_cartels, 1:4, :An)
        Bn, vdA = analyse_links(new_adj, new_adj_cartels, 1:4, :Bn)

        Add = sort(map(sum, eachrow(new_adj)))
        Bdd = sort(map(sum, eachrow(new_adj)))
        return ham, nc, An, Bn, vdham, vdnc, vdA, vdB, Add, Bdd
    end

    hamA = zeros(2, 2)
    ncA = zeros(2, 2)
    AnA = zeros(2, 2)
    BnA = zeros(2, 2)
    nocartel_dd = zeros(10)
    cartel_dd = zeros(10)
    n = 1
    while n <= 50
        @show n
        ham, nc, An, Bn, vdham, vdnc, vdA, vdB, Add, Bdd = apply_method()
        if all((vdham, vdnc, vdA, vdB))
            hamA .= hamA .+ (ham .- hamA) ./ n
            ncA .= ncA .+ (nc .- ncA) ./ n
            AnA .= AnA .+ (An .- AnA) ./ n
            BnA .= BnA .+ (Bn .- BnA) ./ n
            nocartel_dd .= nocartel_dd .+ (Add .- nocartel_dd) ./ n
            cartel_dd .= cartel_dd .+ (Bdd .- cartel_dd) ./ n
            n += 1
        end
    end
    return hamA, ncA, AnA, BnA, nocartel_dd, cartel_dd
end

hamA, ncA, AnA, BnA, nocartel_dd, cartel_dd = main()


function plot_embedded(g, embedding)
    p = plot(proj=:polar, yticks=Int64[])
    AM = adjacency_matrix(g)
    for ij in CartesianIndices(AM)
        i, j = Tuple(ij)
        if AM[i, j] > 0
            esrc = embedding[i]
            edest = embedding[j]
            plot!(p, [esrc[2], edest[2]], [esrc[1], edest[1]], colour=:black, label="", arrow=true)
        end
    end
    # for (src, dests) in enumerate(g.fadjlist)
    #     esrc = embedding[src]
    #     for dest in dests
    #         edest = embedding[dest]
    #         plot!(p, [esrc[2], edest[2]], [esrc[1], edest[1]], colour=:black, label="")
    #     end
    # end
    scatter!(reverse.(embedding), label="", series_annotations=text.(1:length(embedding), :bottom))
    return p
end
