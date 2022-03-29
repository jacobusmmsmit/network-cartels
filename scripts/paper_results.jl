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
        new_adj[player, :] .= 1
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

function analyse_navigability(A, B, cartel, dm)
    nocartel = setdiff(1:size(A, 1), cartel)
    gA = SimpleDiGraph(A)
    gB = SimpleDiGraph(B)
    navdiff(sub1, sub2) = navigability(gB, dm, 10000, sub1, sub2)[1] - navigability(gA, dm, 10000, sub1, sub2)[1]
    cc = navdiff(cartel, cartel)
    cnc = navdiff(cartel, nocartel)
    ncc = navdiff(nocartel, cartel)
    ncnc = navdiff(nocartel, nocartel)
    res = [cc cnc
        ncc ncnc]
    return res
end

function main()
    Random.seed!(2)

    function apply_method()
        N = 10
        toy, (toy_embedding, _) = generate_hyperbolic_graph(N, 0.7, 1.0, 0.15)
        cco = global_clustering_coefficient(toy)
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
            suspiciousness[i] = sum(old_strat .< new_strat)
        end
        suspiciousness .= suspiciousness ./ sum(suspiciousness)
        ham, vdham = analyse_links(new_adj, new_adj_cartels, 1:4, :hamming)
        nc, vdnc = analyse_links(new_adj, new_adj_cartels, 1:4, :nconnections)
        An, vdB = analyse_links(new_adj, new_adj_cartels, 1:4, :An)
        Bn, vdA = analyse_links(new_adj, new_adj_cartels, 1:4, :Bn)
        navdiff = analyse_navigability(new_adj, new_adj_cartels, 1:4, embedded_distance_matrix)
        Odd = sort(map(sum, eachrow(adjacency_matrix(toy))), rev=true)
        Add = sort(map(sum, eachrow(new_adj)), rev=true)
        Bdd = sort(map(sum, eachrow(new_adj_cartels)), rev=true)

        ccnc = global_clustering_coefficient(SimpleDiGraph(new_adj))
        ccc = global_clustering_coefficient(SimpleDiGraph(new_adj_cartels))
        return ham, nc, An, Bn, vdham, vdnc, vdA, vdB, Odd, Add, Bdd, cco, ccnc, ccc, toy, new_adj, toy_embedding, suspiciousness, navdiff
    end

    # Collect data
    hamA = zeros(2, 2)
    ncA = zeros(2, 2)
    AnA = zeros(2, 2)
    BnA = zeros(2, 2)
    navigability_difference = zeros(2, 2)
    clustering_coeff_original = 0.0
    clustering_coeff_nocartel = 0.0
    clustering_coeff_cartel = 0.0
    original_dd = Int64[]
    nocartel_dd = Int64[]
    cartel_dd = Int64[]
    n_cartel_identified = Int64[]
    null_identified = Int64[]
    n = 1
    while n <= 150
        @show n
        ham, nc, An, Bn, vdham, vdnc, vdA, vdB, Odd, Add, Bdd, cco, ccnc, ccc, _, _, _, suspiciousness, navdiff = apply_method()
        if all((vdham, vdnc, vdA, vdB))
            # Calculate rolling average
            hamA .= hamA .+ (ham .- hamA) ./ n
            ncA .= ncA .+ (nc .- ncA) ./ n
            AnA .= AnA .+ (An .- AnA) ./ n
            BnA .= BnA .+ (Bn .- BnA) ./ n
            clustering_coeff_original = clustering_coeff_original + (cco - clustering_coeff_original) / n
            clustering_coeff_nocartel = clustering_coeff_nocartel + (ccnc - clustering_coeff_nocartel) / n
            clustering_coeff_cartel = clustering_coeff_cartel + (ccc - clustering_coeff_cartel) / n
            navigability_difference = navigability_difference + (navdiff - navigability_difference) / n

            append!(original_dd, Odd)
            append!(nocartel_dd, Add)
            append!(cartel_dd, Bdd)
            println(round.(suspiciousness, sigdigits=2))
            append!(n_cartel_identified, sum(1:4 .∈ sortperm(suspiciousness, rev=true)[1:4]))
            append!(null_identified, sum(1:4 .∈ shuffle(1:10)[1:4]))
            n += 1
        end
    end
    ham, nc, An, Bn, vdham, vdnc, vdA, vdB, Odd, Add, Bdd, cco, ccnc, ccc, original_graph, BR_graph, embedding, _, _ = apply_method()
    return hamA, ncA, AnA, BnA, original_dd, nocartel_dd, cartel_dd, clustering_coeff_original, clustering_coeff_nocartel, clustering_coeff_cartel, original_graph, SimpleDiGraph(BR_graph), embedding, n_cartel_identified, null_identified, navigability_difference
end

# Run main
hamA, ncA, AnA, BnA, original_dd, nocartel_dd, cartel_dd, clustering_coeff_original, clustering_coeff_nocartel, clustering_coeff_cartel, original_graph, BR_graph, embedding, n_cartel_identified, null_identified, navigability_difference = main()

mean(n_cartel_identified)
mean(null_identified)

# Table 1: below
navigability_difference

## Analyse data:
# Plot before and after

# Structural changes:
clustering_coeff_original
clustering_coeff_nocartel
clustering_coeff_cartel

function plot_embedded(g, embedding, col="#ffa600", cartel=Int64[], cartel_colour=:red)
    p = plot(proj=:polar, yticks=Int64[], axis=false, gridalpha=0.5)
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
    if length(cartel) == 0
        scatter!(reverse.(embedding), label="", colour=col)
    else
        not_cartel_ppl = embedding[setdiff(1:length(embedding), cartel)]
        cartel_ppl = embedding[cartel]
        scatter!(reverse.(not_cartel_ppl), label="", colour=col)
        scatter!(reverse.(cartel_ppl), label="", colour=cartel_colour)
    end
    return p
end


begin
    # Figure 2
    new_adj_cartels = best_response_dynamics(adjacency_matrix(original_graph), shortestpath_distance_matrix, embedded_distance_matrix, edge_cost, n=2, cartels=[[1, 2, 3, 4]])
    cartel_graph = SimpleDiGraph(new_adj_cartels)
    plot_embedded(cartel_graph, embedding)
    # navigability(original_graph, get_dm(embedding, disth))
    # navigability(BR_graph, get_dm(embedding, disth))
    # navigability(cartel_graph, get_dm(embedding, disth))
    p1 = plot_embedded(original_graph, embedding, "#ffa600")
    plot!(title="Original", titlefontsize=12)
    p2 = plot_embedded(BR_graph, embedding, "#58508d")
    plot!(title="Strategic", titlefontsize=12)
    p3 = plot_embedded(cartel_graph, embedding, "#96db00", collect(1:4))
    plot!(title="Cartels", titlefontsize=12)
    graph_plots = plot(p1, p2, p3, size=(250, 750), layout=(3, 1))
    ylms = (minimum((minimum(original_dd), minimum(nocartel_dd), minimum(cartel_dd))), maximum((maximum(original_dd), maximum(nocartel_dd), maximum(cartel_dd))) + 1)
    podd = scatter(sort(original_dd, rev=true), ylims=ylms)
    plot!(title="⟨k⟩ = $(round(mean(original_dd), sigdigits = 2))", titlefontsize=12)
    pncdd = scatter(sort(nocartel_dd, rev=true), ylims=ylms)
    plot!(title="⟨k⟩ = $(round(mean(nocartel_dd), sigdigits = 2))", titlefontsize=12)
    pcdd = scatter(sort(cartel_dd, rev=true), ylims=ylms)
    plot!(title="⟨k⟩ = $(round(mean(cartel_dd), sigdigits = 2))", titlefontsize=12)
    degree_plots = plot(podd, pncdd, pcdd, layout=(3, 1), size=(250, 750), legend=:none)
    plot(graph_plots, degree_plots, layout=(1, 2), size=(500, 750))
    savefig("out/before_after.pdf")
end

function greedy_navigate(source, destination, graph, distm; verbose=true)
    success = true
    cur = source
    prev = -1
    route = [cur]
    i = 1
    while cur != destination && i <= size(distm, 1)
        i += 1
        possible_next = setdiff(graph.fadjlist[cur], cur)
        if length(possible_next) == 0
            if verbose
                print("Cannot navigate from $source to $destination: ")
                println("$cur is an absorbing node.")
            end
            success = false
            break
        end
        next = possible_next[findmin(pn -> distm[pn, destination], possible_next)[2]]
        if next == prev
            if verbose
                print("Cannot navigate from $source to $destination: ")
                println("$cur and $next form a greedy loop.")
            end
            success = false
            break
        end
        append!(route, next)
        prev, cur = cur, next
    end
    return route, success
end

function navigability(digraph, dm, N=100000, from_subset=1:size(dm, 1), to_subset=1:size(dm, 1))
    success_rate = 0.0
    mean_length = 0.0
    for i in 1:N
        s = sample(from_subset)
        t = sample(setdiff(to_subset, s))
        route, success = greedy_navigate(s, t, digraph, dm, verbose=false)
        lr = length(route)
        success_val = Int(success)
        success_rate = success_rate + (success_val - success_rate) / i
        mean_length = mean_length + (lr - mean_length) / i
    end
    mean_length, success_rate
end

