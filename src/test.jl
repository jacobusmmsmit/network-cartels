using Graphs
using SimpleWeightedGraphs
# using GraphIO
using StatsBase
using Plots
using IterTools
using Distributions
using Measures

# p2p = loadgraph("data/p2p-Gnutella31.txt", "graph_key", EdgeListFormat())

include("game_theory.jl")
include("generate_HRG.jl")
include("helper_functions.jl")
include("routing.jl")


function ns(npeers=4, α=0.2)
    # peers = [rand(Uniform(0, 100), 2) for _ in 1:npeers]
    # sort!(peers, by=first)
    # dm = [diste(p1, p2) for (p1, p2) in Iterators.product(peers, peers)]
    # strategies = [sample(1:npeers, rand(1:k), replace=false) for k in 1:npeers]
    # g = construct_graph(strategies, dm)
    g, (peers, edges) = generate_hyperbolic_graph(npeers, 0.7, 1.0, 0.2)
    dm = [diste(p1, p2) for (p1, p2) in Iterators.product(peers, peers)]
    ds = floyd_warshall_shortest_paths(g)
    sps = ds.dists # from i to j
    @show edges
    strategies = VE_to_strategies(peers, edges)

    function cfn(player, strategies, distm,
        g=construct_graph(strategies, distm),
        ds=floyd_warshall_shortest_paths(g),
        sps=ds.dists) # from i to j
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
        p = plot(legend=false, aspect_ratio=1.0)
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

    # begin
    #     p = plot(legend=false, aspect_ratio=1.0, proj=:polar, ticks=false, ylabel="Foo", xlabel="Bar", margin=3mm)
    #     for (s, t) in edges
    #         @show (s, t)
    #         plot!([s[2], t[2]], [s[1], t[1]], lc=:gray, lw=0.5)
    #     end
    #     p
    #     scatter!(last.(peers), first.(peers), ms=3.0) # series_annotations=text.(1:length(peers), :bottom)
    # end

    new_strategies = best_response_dynamics(cfn, strategies, dm; n=1)
    NE_costs = evaluate_cost(cfn, 1:npeers, new_strategies, dm)
    NE_social_cost = sum(NE_costs)
    @show NE_costs
    @show NE_social_cost
    begin
        p2 = plot(legend=false, aspect_ratio=1.0)
        scatter!(first.(peers), last.(peers), series_annotations=text.(1:length(peers), :bottom))
        wts = mapreduce(vcat, enumerate(new_strategies)) do x
            i = x[1]
            ds = x[2]
            map(dest -> dm[i, dest], ds)
        end
        sources = inverse_rle(1:npeers, length.(new_strategies))
        dests = collect(Iterators.flatten(new_strategies))
        for (_s, _d, w) in zip(sources, dests, wts)
            s = peers[_s]
            d = peers[_d]
            midpoint = (s[1] + d[1]) / 2, (s[2] + d[2]) / 2
            plot!([s[1], midpoint[1]], [s[2], midpoint[2]], colour=:black, arrow=true, label="")
            plot!([midpoint[1], d[1]], [midpoint[2], d[2]], colour=:black, label="")
        end
        p2
    end

    nav_score = navigability(construct_graph(new_strategies, dm), dm, 1000)
    # return p
    @show nav_score

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
    plot(p, p2, size=(1000, 500))
end

# main(5, 1.0)

ns(4)
function VE_to_strategies(V, E)
    VDict = Dict(V .=> 1:length(V))
    dsts = get.(Ref(VDict), last.(E), 0)
    srcs = get.(Ref(VDict), first.(E), 0)
    strategies = [Int64[] for _ in 1:length(V)]
    for (s, t) in zip(srcs, dsts)
        push!(strategies[s], t)
    end
    return strategies
end



# function main()
#     N = 9
#     ms = zeros(N)
#     ml = zeros(N)
#     for i in 2:N
#         @show i
#         success_rate = 0.0
#         mean_length = 0.0
#         for n in 1:100
#             lr, suc = ns(i, 1.0)
#             success_rate = success_rate + (suc - success_rate)/n
#             mean_length = mean_length + (lr - mean_length)/n
#         end
#         ml[i] = mean_length
#         ms[i] = success_rate
#     end
#     ms, ml
# end
# ms, ml = main()
# plot(ml)

## TODO: Implement check for Nash Equilibrium
## TODO: Implement sample peers from hyperbolic space
## TODO: Implement hyperbolic distances
## TODO: Decide on inputs such as
# number of best response dynamics,
# number of peers,
# space size,
# null model
