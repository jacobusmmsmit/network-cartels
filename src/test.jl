using Graphs
using SimpleWeightedGraphs
using GraphIO
using StatsBase
using Plots

# p2p = loadgraph("data/p2p-Gnutella31.txt", "graph_key", EdgeListFormat())

norme(x) = √((x[1])^2 + (x[2])^2)
diste(x, y) = norme(x .- y)
disth(x, y) = acosh(1 + (2*diste(x, y)^2)/((1 - norme(x)^2)*(1 - norme(y)^2)))

function construct_graph(strategies, distm)
    sources = inverse_rle(1:5, length.(strategies))
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

function main(npeers = 5, α = 1.0)
    peers = [rand(2) for _ in 1:npeers]
    dm = [diste(p1, p2) for (p1, p2) in Iterators.product(peers, peers)]
    fn = (player, edges, sps) -> p2p_cost(player, edges, sps, dm, α)

    strategies = [sample(1:5, rand(1:5), replace = false) for _ in 1:npeers]
    g = construct_graph(strategies, dm)
    ds = floyd_warshall_shortest_paths(g)
    sps = ds.dists # from i to j

    costs = evaluate_cost(fn, 1:5, strategies, sps)

    scatter(first.(peers), last.(peers))
    wts = mapreduce(vcat, enumerate(strategies)) do x
        i = x[1]
        ds = x[2]
        map(dest -> dm[i, dest], ds)
    end
    sources = inverse_rle(1:5, length.(strategies))
    dests = collect(Iterators.flatten(strategies))
    @show length.((sources, dests, wts))
    for (s, d, w) in zip(sources, dests, wts)
        peer_s = peers[]
        @show s
        @show d
        @show w
        plot!([s[1], d[1]], [s[2], d[2]])
    end
    plot
end

main()
