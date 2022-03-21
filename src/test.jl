using Graphs
using SimpleWeightedGraphs
using GraphIO
using StatsBase
using Plots

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

function evaluate_cost(fn, players, strategies, shortest_paths)
    map(zip(players, strategies)) do (player, dests)
        edges = [(player, dest) for dest in dests]
        fn(player, edges, shortest_paths)
    end
end

function p2p_cost(player, strategy, sps, distm, α)
    cost = α * length(strategy)
    for other in Iterators.flatten((1:player-1, (player+1):(size(distm, 1))))
        cost += sps[player, other] / distm[player, other]
    end
    cost
end

# UnitRangeWithOut
ur_wo(ur, without) = Iterators.flatten((ur.start:(without-1), (without+1):ur.stop))

function greedy_navigate(source, destination, graph, distm)
    cur = source
    prev = -1
    route = [cur]
    while cur != destination
        possible_next = setdiff(graph.weights[:, cur].nzind, cur)
        if length(possible_next) == 0
            break
        end
        # println(map(pn -> distm[pn, destination], possible_next))
        # @show possible_next
        # @show distm
        next = possible_next[findmin(pn -> distm[pn, destination], possible_next)[2]]
        if next == prev
            @show next
            @show prev
            break
        end
        append!(route, next)
        prev, cur = cur, next
    end
    route
end

function main(npeers=4, α=0.2)
    peers = [rand(2) for _ in 1:npeers]
    sort!(peers, by=first)
    dm = [diste(p1, p2) for (p1, p2) in Iterators.product(peers, peers)]
    fn = (player, edges, sps) -> p2p_cost(player, edges, sps, dm, α)

    strategies = [sample(1:npeers, rand(1:npeers), replace=false) for _ in 1:npeers]
    g = construct_graph(strategies, dm)
    ds = floyd_warshall_shortest_paths(g)
    sps = ds.dists # from i to j

    costs = evaluate_cost(fn, 1:npeers, strategies, sps)
    social_cost = sum(costs)

    println(greedy_navigate(1, 2, g, dm))

    begin
        p = plot(legend=false)
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
end

main(10)
