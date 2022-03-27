export construct_graph, evaluate_cost, p2p_cost, best_response, best_response_dynamics, social_optimum

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