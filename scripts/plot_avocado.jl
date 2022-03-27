include("../src/NetworkCartelDetection.jl")
using .NetworkCartelDetection
using Plots
using Measures
using ColorSchemes
using Random
gr()

Random.seed!(6)

function main(; N=50, α=0.8, ν=1.0, T=0.3)
    g, (peers, edges) = generate_hyperbolic_graph(N, α, ν, T)
    p = plot(aspect_ratio=1.0, proj=:polar, ticks=false, ylabel="Foo", xlabel="Bar", margin=3mm, legend=:outertop)
    R = 2log(N / ν)
    @show R
    # # Testing purposes:
    # ρ = LinRange(0.0, R, 50)
    # θ = LinRange(0.0, 2π, 360)
    # Actual values:
    ρ = LinRange(0.0, R, 500)
    θ = LinRange(0.0, 2π, 3600)

    pt = peers[findmin(x -> abs(x[1] - (2R / 3)), peers)[2]] # Centre of ball
    # pt = (R / 2, 2)
    # contourf currently doesn't work with GR backend, so we have to cheat...
    funp(ρ, θ) = disth(pt, (ρ, θ)) < R ? (isapprox(disth(pt, (ρ, θ)), R, atol=2e-1) ? 0 : floor(disth(pt, (ρ, θ)))) : 1.1R
    # funp(ρ, θ) = floor(disth(pt, (ρ, θ)))
    heatmap!(θ, ρ, (θ, ρ) -> funp(ρ, θ), proj=:polar, colorbar=false, color=cgrad(:tempo, rev=true))
    # plot!([2, 2 + π], [R / 2, R / 2], lw=2, lc=:black)
    # plot!(title = "Distance from Point Half-way to Radius") # put the title as caption
    first_special = true
    first_normal = true
    for (s, t) in edges
        if s == pt || t == pt
            if first_special
                first_special = false
                plot!([s[2], t[2]], [s[1], t[1]], lc=:red, lw=1.0, label="Ball centre's edges")
            else
                plot!([s[2], t[2]], [s[1], t[1]], lc=:red, lw=1.0, label="")
            end
        else
            if first_normal
                first_normal = false
                plot!([s[2], t[2]], [s[1], t[1]], lc=:gray, lw=0.5, label="Edges")
            else
                plot!([s[2], t[2]], [s[1], t[1]], lc=:gray, lw=0.5, label="")
            end

        end
    end
    scatter!(last.(peers), first.(peers), mc=:blue, ms=3.0, label="Nodes") # series_annotations=text.(1:length(peers), :bottom)
    scatter!([0], [0], mc=:black, label="Origin", shape=:plus, ms=5, msw=3.0)
    scatter!([pt[2]], [pt[1]], mc=:white, label="Ball Centre", shape=:diamond, msw=2.0)
    plot!(size=(300, 300), yticks=[0], legend=:none) # set yticks instead of hiding as this bugs the plot

    # The coolest thing I've made by accident:
    # let
    #     N = 100
    #     ν = 1.0
    #     R = 2log(N/ν)
    #     ρ = LinRange(0., 7, 200)
    #     θ = LinRange(0., 2π, 360)
    #     funp(ρ, θ) = round(disth((3, 0), (ρ, θ)))
    #     Plots.heatmap(ρ, θ, funp, proj = :polar)
    # end
    savefig("out/hyperbolic_ball_with_nodes_nolegend.pdf")
    return p
end

main(N=100, α=0.7, ν=0.5, T=0.0)