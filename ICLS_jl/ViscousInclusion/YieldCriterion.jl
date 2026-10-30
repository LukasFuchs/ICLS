# ======================================================================= #
# YieldCriterion.jl
#
# Onset of yielding for an isolated elliptical inclusion (inherited weak zone,
# IWZ) in an isoviscous matrix, using the closed-form amplification factors
# and a Drucker-Prager yield stress with lithostatic pressure (as in LaMEM
# with p_litho_plast = 1). No numerical model is solved.
#
# Quantities
#   Δη  = η_IWZ / η_m                    viscosity contrast
#   Λ   = τ_y,IWZ / τ_y,m                strength ratio (same depth)
#   τ_bg* / τ_y,m                        background stress at which yielding
#                                        starts, relative to the matrix strength
#   inclusion yields first : τ_bg*/τ_y,m = Λ / A_τ
#   matrix at the interface: τ_bg*/τ_y,m = 1 / max(A_ε, A_τ)
#   regime boundary        : Λ = Δη (Δη ≤ 1),  Λ = 1 (Δη ≥ 1)
#
# Figure (output/Yield/YieldCriterion.png):
#   (a) τ_bg*/τ_y,m over (Δη, Λ) for one shape and orientation
#   (b) τ_bg*/τ_y,m versus β for several a/b (fixed Δη, Λ)
#   (c) Λ versus depth for the IWZ parameters of the LaMEM multi-IWZ setup
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. ViscousInclusion/YieldCriterion.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "InclusionTools.jl"))
using Plots, Printf

function run_yield()
    IT = InclusionTools
    # ===================== Parameters ================================== #
    t_a     =   3.0                         # (a) aspect ratio a/b
    β_a     =   45.0                        # (a) angle long axis - shortening direction [deg]
    Δη_b    =   0.1                         # (b) viscosity contrast
    Λ_b     =   0.6                         # (b) strength ratio
    ratios  =   [1.0, 2.0, 3.0, 5.0, 8.0]   # (b) aspect ratios
    # (c) matrix and IWZ parameters (LaMEM multi-IWZ setup, crustal phases)
    c_m, φ_m =  20e6, 15.0                  # matrix cohesion [Pa], friction angle [deg]
    iwz     =   [(8.03e6, 8.67591), (6.75e6, 11.89), (4.60e6, 7.90707), (5.94e6, 11.4435),
                 (9.68e6, 12.1261), (6.61e6, 8.64361), (10.2e6, 11.4755)]
    ρ, g    =   2700.0, 9.81                # crustal density [kg/m³], gravity [m/s²]
    zmax    =   40e3                        # [m]
    outdir  =   joinpath(@__DIR__, "output", "Yield")
    # =================================================================== #
    mkpath(outdir)
    cmap    =   cgrad(:viridis)
    fs      =   (framestyle = :box, titlefontsize = 13, guidefontsize = 13, tickfontsize = 11)
    # ------------------------------------------------------------------- #
    # (a) onset over (Δη, Λ), plotted in log10 coordinates
    # ------------------------------------------------------------------- #
    lΔη     =   collect(-3.0:0.0125:3.0)
    lΛ      =   collect(-1.0:0.00625:1.0)
    R       =   [IT.yield_onset(10.0^x, t_a, β_a, 10.0^y).ratio for y in lΛ, x in lΔη]
    Z       =   log10.(R)
    zlim    =   (log10(0.03), 0.0)
    lev     =   [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    pa  =   heatmap(lΔη, lΛ, Z; color = cmap, clims = zlim, colorbar = false,
                    xlims = extrema(lΔη), ylims = extrema(lΛ),
                    xticks = (-3:3, ["10⁻³", "10⁻²", "10⁻¹", "1", "10¹", "10²", "10³"]),
                    yticks = (-1:1, ["10⁻¹", "1", "10¹"]),
                    xlabel = "\$\\Delta\\eta = \\eta_{IWZ}/\\eta_m\$",
                    ylabel = "\$\\Lambda = \\tau_{y,IWZ}/\\tau_{y,m}\$",
                    title = "(a) far-field stress at first yielding\n(a/b = $(t_a), β = $(β_a)°)",
                    legend = false, fs...)
    contour!(pa, lΔη, lΛ, Z; levels = log10.(lev), color = :black, lw = 0.6, colorbar_entry = false)
    # contour labels: where the contours cross a row near the top (vertical
    # contours, Λ > Δη) and a column near the right edge (Λ < 1, stiff side)
    irow    =   findfirst(>=(0.93), lΛ)
    jcol    =   findfirst(>=(2.6), lΔη)
    for v in lev
        lv  =   log10(v)
        j   =   findfirst(j -> (Z[irow, j] - lv) * (Z[irow, j+1] - lv) <= 0, 1:length(lΔη)-1)
        j === nothing || annotate!(pa, lΔη[j], 0.965, text(string(v), 9, :black, :center))
        i   =   findfirst(i -> (Z[i, jcol] - lv) * (Z[i+1, jcol] - lv) <= 0, 1:length(lΛ)-1)
        i === nothing || lΛ[i] > -0.05 || annotate!(pa, 2.75, lΛ[i] + 0.035, text(string(v), 9, :black, :center))
    end
    vline!(pa, [0.0]; color = :gray30, lw = 0.8, ls = :dot)
    hline!(pa, [0.0]; color = :gray30, lw = 0.8, ls = :dot)
    plot!(pa, lΔη, lΔη; color = :black, lw = 1.4, ls = :dash)                  # Λ = Δη
    plot!(pa, lΔη, min.(lΔη, 0.0); color = :red, lw = 2.2)                     # regime boundary
    annotate!(pa, -2.95, 0.80, text("left of dashed line:\nweaker in viscosity than in strength", 10, :left, :black))
    annotate!(pa, -2.95, 0.30, text("red line: which yields first\nabove: matrix at the edges of the zone\nbelow: the zone itself", 10, :left, :black))
    annotate!(pa,  2.95, 0.40, text("right of dashed line:\nweaker in strength\nthan in viscosity", 10, :right, :black))
    annotate!(pa,  1.30, 0.86, text("Λ = Δη", 10, :left, :black))
    # colorbar as its own panel with the actual ratio as labels
    tk      =   [0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
    vv      =   collect(range(zlim..., length = 256))
    pcb =   heatmap([0.0, 1.0], vv, repeat(vv, 1, 2); color = cmap, clims = zlim, colorbar = false,
                    xticks = false, yticks = (log10.(tk), string.(tk)), ymirror = true, xlims = (0, 1),
                    ylims = zlim, framestyle = :box, tickfontsize = 11, titlefontsize = 12,
                    title = "\$\\tau_{bg}^*/\\tau_{y,m}\$", left_margin = 0Plots.mm, right_margin = 8Plots.mm)
    # ------------------------------------------------------------------- #
    # (b) onset versus orientation for several aspect ratios
    # ------------------------------------------------------------------- #
    βs      =   collect(0.0:0.25:90.0)
    regime  =   IT.yield_onset(Δη_b, 3.0, 45.0, Λ_b).regime == :zone ? "zone yields first" : "edges yield first"
    pb  =   plot(; xlims = (0, 90), xticks = 0:15:90, ylims = (0, 1.05), xlabel = "β [°]",
                 ylabel = "\$\\tau_{bg}^*/\\tau_{y,m}\$", legend = :bottomleft, legendfontsize = 10,
                 left_margin = 22Plots.mm,
                 title = "(b) versus orientation\n(Δη = $(Δη_b), Λ = $(Λ_b): $(regime))", fs...)
    hline!(pb, [1.0]; color = :gray, lw = 0.6, label = "")
    colors  =   [:black; palette(:tab10)[2:length(ratios)]]
    for (m, t) in enumerate(ratios)
        plot!(pb, βs, [IT.yield_onset(Δη_b, t, β, Λ_b).ratio for β in βs]; color = colors[m], lw = 1.8,
              label = t == 1.0 ? "circle" : "a/b = $(t)")
    end
    # ------------------------------------------------------------------- #
    # (c) strength ratio versus depth for the IWZ parameters
    # ------------------------------------------------------------------- #
    zs      =   collect(range(0.0, zmax, length = 401))
    pc  =   plot(; yflip = true, xlabel = "Λ", ylabel = "depth [km]", legend = :bottomleft,
                 legendfontsize = 9, left_margin = 6Plots.mm,
                 title = @sprintf("(c) Λ versus depth\n(matrix: c = %g MPa, φ = %g°; ρ = %g kg/m³)", c_m / 1e6, φ_m, ρ),
                 fs...)
    for (k, (c, φ)) in enumerate(iwz)
        Λz  =   [IT.strength_ratio(c, φ, c_m, φ_m, IT.lithostatic_pressure(ρ, g, z)) for z in zs]
        plot!(pc, Λz, zs ./ 1e3; color = palette(:tab10)[k], lw = 1.6,
              label = @sprintf("#%d: c = %.3g MPa, φ = %.1f°", k, c / 1e6, φ))
    end
    # ------------------------------------------------------------------- #
    p   =   plot(pa, pcb, pb, pc; layout = grid(1, 4, widths = [0.37, 0.018, 0.33, 0.282]),
                 size = (2200, 700), left_margin = 10Plots.mm, right_margin = 10Plots.mm,
                 bottom_margin = 10Plots.mm, top_margin = 6Plots.mm)
    savefig(p, joinpath(outdir, "YieldCriterion.png"))

    # console summary
    @printf("\nOnset of yielding, a/b = %g, β = %g°:\n", t_a, β_a)
    @printf("%8s %6s | %12s %8s\n", "Δη", "Λ", "τ_bg*/τ_y,m", "regime")
    for (Δη, Λ) in ((0.01, 0.6), (0.1, 0.6), (0.5, 0.6), (1.0, 0.6), (0.5, 0.3), (10.0, 0.6), (10.0, 2.0))
        o   =   IT.yield_onset(Δη, t_a, β_a, Λ)
        @printf("%8g %6g | %12.3f %8s\n", Δη, Λ, o.ratio, o.regime)
    end
    println("\nFigure saved to $(outdir)")
    return nothing
end

run_yield()
