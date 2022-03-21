using Graphs
using SimpleWeightedGraphs
using GraphIO
using StatsBase
using Plots
using IterTools
using Distributions

# p2p = loadgraph("data/p2p-Gnutella31.txt", "graph_key", EdgeListFormat())

norme(x) = √((x[1])^2 + (x[2])^2)
diste(x, y) = norme(x .- y)
disth(x, y) = acosh(1 + (2 * diste(x, y)^2) / ((1 - norme(x)^2) * (1 - norme(y)^2)))

function construct_graph(strategies, distm)
    sources = inverse_rle(1:length(strategies), length.(strategies))
    destinations = collect(Iterators.flatten(strategies))
    wts = mapreduce(vcat, enumerate(strategies)) do x
        i = x[1]
        i_dests = x[2]
        map(dest -> distm[i, dest], i_dests)
    end
    return SimpleWeightedDiGraph(sources, destinations, wts)
end

function evaluate_cost(cfn, players, strategies, distm)
    g = construct_graph(strategies, distm)
    ds = floyd_warshall_shortest_paths(g)
    sps = ds.dists # from i to j
    return [cfn(player, strategies, distm, g, ds, sps) for player in players]
end

function p2p_cost(player, strategy, sps, distm, α)
    cost = α * length(strategy)
    for other in Iterators.flatten((1:player-1, (player+1):(size(distm, 1))))
        cost += sps[player, other] / distm[player, other]
    end
    return cost
end

function best_response(cost_function, player, strategies, distm)
    players = 1:length(strategies)
    min_cost = evaluate_cost(cost_function, players, strategies, distm)[player]
    BR = strategies[player]
    for strategy in subsets(collect(ur_wo(1:length(strategies), player)))
        # @show strategy
        new_strategies = copy(strategies)
        new_strategies[player] = strategy
        cur_cost = evaluate_cost(cost_function, players, new_strategies, distm)[player]
        if cur_cost < min_cost
            min_cost = cur_cost
            BR = strategy
        end
    end
    return BR, min_cost
end

function best_response_dynamics(cost_function, strategies, distm; n = 1)
    players = 1:length(strategies)
    new_strategies = copy(strategies)
    for _ in 1:n
        for player in players
            new_strategies[player] = best_response(cost_function, player, new_strategies, distm)[1]
        end
    end
    return new_strategies
end

function social_optimum(cost_function, distm, hardcore = false)
    players = 1:size(distm, 1)
    min_cost = Inf
    social_optimum_solution = [Int64[] for _ in players]
    strategies = [Int64[] for _ in players]
    all_possible_strategies = product(
        (collect(subsets(players)) for _ in players)...
        )
    if length(all_possible_strategies) > 65536
        println(length(all_possible_strategies))
        throw(ArgumentError("all possible strategies is too large"))
    end
    for _strategies in all_possible_strategies
        if !all(length.(_strategies) .> 0)
            continue
        end
        strategies .= _strategies
        cur_cost = sum(evaluate_cost(cost_function, players, strategies, distm))
        if cur_cost < min_cost
            social_optimum_solution = strategies
            min_cost = cur_cost
        end
    end
    return social_optimum_solution, min_cost
end

# UnitRangeWithOut
ur_wo(ur, without) = Iterators.flatten((ur.start:(without-1), (without+1):ur.stop))

function greedy_navigate(source, destination, graph, distm; verbose = true)
    success = true
    cur = source
    prev = -1
    route = [cur]
    i = 1
    while cur != destination && i <= size(distm, 1)
        i += 1
        possible_next = setdiff(graph.weights[:, cur].nzind, cur)
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

function navigability(digraph, dm, N = 1000)
    n = size(digraph.weights, 1)
    success_rate = 0.0
    mean_length = 0.0
    for i in 1:N
        s, t = sample(1:n, 2, replace = false)
        route, success = greedy_navigate(s, t, digraph, dm, verbose = false)
        lr = length(route)
        success_val = Int(success)
        success_rate = success_rate + (success_val - success_rate)/i
        mean_length = mean_length + (lr - mean_length)/i
    end
    mean_length, success_rate
end



function ns(npeers=4, α=0.2)
    peers = [rand(Uniform(0, 100), 2) for _ in 1:npeers]
    sort!(peers, by=first)
    dm = [diste(p1, p2) for (p1, p2) in Iterators.product(peers, peers)]
    strategies = [sample(1:npeers, rand(1:k), replace=false) for k in 1:npeers]
    g = construct_graph(strategies, dm)
    ds = floyd_warshall_shortest_paths(g)
    sps = ds.dists # from i to j
    
    function cfn(player, strategies, distm,
            g = construct_graph(strategies, distm),
            ds = floyd_warshall_shortest_paths(g),
            sps = ds.dists) # from i to j
        edges = [(player, dest) for dest in strategies[player]]
        p2p_cost(player, edges, sps, distm, α)
    end

    # fn = (player, edges) -> p2p_cost(player, edges, sps, dm, α)
    costs = evaluate_cost(cfn, 1:npeers, strategies, dm)
    social_cost = sum(costs)
    # @show costs
    # @show social_cost

    # println(greedy_navigate(1, 2, g, dm))

    begin
        p = plot(legend=false, aspect_ratio = 1.0)
        scatter!(first.(peers), last.(peers), series_annotations=text.(1:length(peers), :bottom))
        wts = mapreduce(vcat, enumerate(strategies)) do x
            i = x[1]
            ds = x[2]
            map(dest -> dm[i, dest], ds)
        end
        sources = inverse_rle(1:npeers, length.(strategies))
        dests = collect(Iterators.flatten(strategies))
        for (_s, _d, w) in zip(sources, dests, wts)
            s = peers[_s]
            d = peers[_d]
            midpoint = (s[1] + d[1]) / 2, (s[2] + d[2]) / 2
            plot!([s[1], midpoint[1]], [s[2], midpoint[2]], colour=:black, arrow=true, label="")
            plot!([midpoint[1], d[1]], [midpoint[2], d[2]], colour=:black, label="")
        end
        p
    end
    new_strategies = best_response_dynamics(cfn, strategies, dm; n = 1)
    NE_costs = evaluate_cost(cfn, 1:npeers, new_strategies, dm)
    NE_social_cost = sum(NE_costs)
    # @show NE_costs
    # @show NE_social_cost
    # begin
    #     p2 = plot(legend=false, aspect_ratio = 1.0)
    #     scatter!(first.(peers), last.(peers), series_annotations=text.(1:length(peers), :bottom))
    #     wts = mapreduce(vcat, enumerate(new_strategies)) do x
    #         i = x[1]
    #         ds = x[2]
    #         map(dest -> dm[i, dest], ds)
    #     end
    #     sources = inverse_rle(1:npeers, length.(new_strategies))
    #     dests = collect(Iterators.flatten(new_strategies))
    #     for (_s, _d, w) in zip(sources, dests, wts)
    #         s = peers[_s]
    #         d = peers[_d]
    #         midpoint = (s[1] + d[1]) / 2, (s[2] + d[2]) / 2
    #         plot!([s[1], midpoint[1]], [s[2], midpoint[2]], colour=:black, arrow=true, label="")
    #         plot!([midpoint[1], d[1]], [midpoint[2], d[2]], colour=:black, label="")
    #     end
    #     p2
    # end

    nav_score = navigability(construct_graph(new_strategies, dm), dm, 1000)
    return nav_score
    # @show nav_score
    
    # soc_opt_strats, soc_opt_cost = social_optimum(cfn, dm)
    # begin
    #     p3 = plot(legend=false, aspect_ratio = 1.0)
    #     scatter!(first.(peers), last.(peers), series_annotations=text.(1:length(peers), :bottom))
    #     wts = mapreduce(vcat, enumerate(soc_opt_strats)) do x
    #         i = x[1]
    #         ds = x[2]
    #         map(dest -> dm[i, dest], ds)
    #     end
    #     sources = inverse_rle(1:npeers, length.(soc_opt_strats))
    #     dests = collect(Iterators.flatten(soc_opt_strats))
    #     for (_s, _d, w) in zip(sources, dests, wts)
    #         s = peers[_s]
    #         d = peers[_d]
    #         midpoint = (s[1] + d[1]) / 2, (s[2] + d[2]) / 2
    #         plot!([s[1], midpoint[1]], [s[2], midpoint[2]], colour=:black, arrow=true, label="")
    #         plot!([midpoint[1], d[1]], [midpoint[2], d[2]], colour=:black, label="")
    #     end
    #     p3
    # end
    # @show soc_opt_cost
    # plot(p, p2, size = (1000, 500))
end

# main(5, 1.0)

function main()
    N = 9
    ms = zeros(N)
    ml = zeros(N)
    for i in 2:N
        @show i
        success_rate = 0.0
        mean_length = 0.0
        for n in 1:100
            lr, suc = ns(i, 1.0)
            success_rate = success_rate + (suc - success_rate)/n
            mean_length = mean_length + (lr - mean_length)/n
        end
        ml[i] = mean_length
        ms[i] = success_rate
    end
    ms, ml
end
ms, ml = main()
plot(ml)

## TODO: Implement check for Nash Equilibrium
## TODO: Implement sample peers from hyperbolic space
## TODO: Implement hyperbolic distances
## TODO: Decide on inputs such as
# number of best response dynamics,
# number of peers,
# space size,
# null model
