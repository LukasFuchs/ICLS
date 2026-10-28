# ======================================================================= #
# ParameterMaps.jl
#
# Analytical parameter maps for an isolated elliptical inclusion
# (Schmid & Podladchikov, 2003): inclusion strain rate ε̇_II, deviatoric
# stress τ_II, and dissipative energy ψ = τ_II ⋅ ε̇_II as a function of the
# orientation θ, the viscosity contrast Δη = ηc/ηm, and the aspect ratio a/b.
# No numerical model is solved; inside the inclusion all fields are
# homogeneous, so the values are exact.
#
# Figures (in output/Maps/):
#   1. Maps_<deformation>.png   : log10 of ε̇_II, τ_II, ψ normalized by the
#                                 background, over (Δη, θ), one column per a/b
#   2. ShapeFactor_<...>.png    : effect of shape and orientation, i.e. the
#                                 value relative to a circular inclusion with
#                                 the same Δη. At fixed Δη this factor is the
#                                 same for ε̇_II and τ_II, and squared for ψ.
#   3. Summary_<...>.png        : weak and rigid limits, maximum of ψ and
#                                 the contrast at which it occurs, vs θ
#   4. Amplification_<...>.png  : amplification factors A_ε = ε̇_II/ε̇_bg,
#                                 A_τ = τ_II/τ_bg, A_ψ = ψ/ψ_bg (circle included),
#                                 contours at fixed values of the factor
#   5. Decomposition.png        : A_ε = C(Δη) · S(Δη, t, β): contrast part
#                                 C = 2/(1+Δη), shape/orientation part S in the
#                                 weak and rigid limits, and A_ε as their product
#                                 (closed form, independent of the deformation type)
#   6. AmpShape_ab_beta.png     : A_ε and S as maps over orientation β and
#                                 aspect ratio a/b, one column per contrast Δη
#                                 (C is constant per column and given in the title;
#                                 A_τ = Δη A_ε and A_ψ = Δη A_ε² have the same pattern)
#   7. AmpBand.png              : A_ε vs Δη for all a/b in one panel; the band
#                                 between β = 0° and 45° contains every orientation,
#                                 since A_ε(β)² = A_ε(0)² cos²2β + A_ε(45°)² sin²2β
#   + CSV files with the maps (= amplification factors) and the summary curves.
#
# Notes
#   - For pure shear or simple shear alone, the inclusion invariants repeat
#     every 90° and are symmetric about 45°; simple shear equals pure shear
#     shifted by 45°. θ = 0 means the long axis is parallel to x.
#   - Background normalization: ε̇_bg, τ_bg = 2 ηm ε̇_bg, ψ_bg = τ_bg ε̇_bg.
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. ViscousInclusion/ParameterMaps.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "InclusionTools.jl"))
using Plots, Printf

function run_maps()
    IT = InclusionTools
    # ===================== Parameters ================================== #
    deformation =   :PureShear                  # :PureShear, :SimpleShear, or (ε̇=..., γ̇=...)
    ratios      =   [2.0, 3.0, 5.0, 8.0]        # aspect ratios a/b (a circle, a/b = 1, is the reference)
    θs          =   collect(0.0:0.5:90.0)       # orientation [deg]
    logΔη       =   collect(-5.0:0.025:5.0)     # log10(ηc/ηm) for the maps
    logΔη_fine  =   collect(-4.0:0.001:4.0)     # finer range for the position of max ψ
    save_csv    =   true
    outdir      =   joinpath(@__DIR__, "output", "Maps")
    # =================================================================== #
    mkpath(outdir)
    ε̇, γ̇       =   IT.flow_rates(deformation, 1.0)     # ε̇_bg = 1 -> τ_bg = 2, ψ_bg = 2
    ref         =   (εII = 1.0, τII = 2.0, ψ = 2.0)
    dname       =   deformation isa Symbol ? string(deformation) : "GeneralShear"

    # normalized inclusion values (ε̇_II, τ_II, ψ) for aspect ratio t, θ, Δη
    function incvals(t, θ, Δη)
        s   =   IT.EllipseSolution(a = t, b = 1.0, θ = θ, ηm = 1.0, ηc = Δη, ε̇ = ε̇, γ̇ = γ̇)
        v   =   IT.inclusion_values(s)
        return v.εII / ref.εII, v.τII / ref.τII, v.ψ / ref.ψ
    end
    maps(t) = begin
        E   =   zeros(length(θs), length(logΔη)); T = similar(E); Ψ = similar(E)
        for (i, θ) in enumerate(θs), (j, l) in enumerate(logΔη)
            E[i, j], T[i, j], Ψ[i, j] = incvals(t, θ, 10.0^l)
        end
        (ε = E, τ = T, ψ = Ψ)
    end

    circle  =   maps(1.0)
    M       =   Dict(t => maps(t) for t in ratios)

    # ------------------------------------------------------------------- #
    # Plot helpers. Each row of maps gets its own colorbar as a separate,
    # narrow panel: GR ignores colorbar_ticks, and a colorbar attached to the
    # last map squeezes that map. Contours are drawn on the same log10 field
    # as the colors, because GR drops contour levels outside the color limits.
    # ------------------------------------------------------------------- #
    cmap    =   cgrad(:RdBu, rev = true)
    xl      =   "\$\\log_{10}\\,\\Delta\\eta\$"
    tlab(v) =   [replace(string(round(x; sigdigits = 3)), r"\.0$" => "") for x in v]
    function cbar(lim, ticks, labels, title)
        v   =   collect(range(-lim, lim, length = 256))
        heatmap([0.0, 1.0], v, repeat(v, 1, 2); color = cmap, clims = (-lim, lim), colorbar = false,
                xticks = false, yticks = (ticks, labels), ymirror = true, xlims = (0, 1),
                ylims = (-lim, lim), framestyle = :box, ylabel = title, guidefontsize = 9,
                tickfontsize = 7, right_margin = 6Plots.mm)
    end
    # nr rows of nc maps, each row followed by its colorbar
    maplayout(nr, nc; wcb = 0.035) = grid(nr, nc + 1, widths = vcat(fill((1 - wcb) / nc, nc), wcb))

    # ------------------------------------------------------------------- #
    # Figure 1: absolute values, diverging colors centered at background
    # ------------------------------------------------------------------- #
    common  =   (framestyle = :box, xlims = extrema(logΔη), ylims = extrema(θs), colorbar = false,
                 yticks = 0:15:90, titlefontsize = 10, guidefontsize = 9, tickfontsize = 7)
    rows    =   ((:ε, "\$\\log_{10}\\,\\dot{\\varepsilon}_{II}/\\dot{\\varepsilon}_{bg}\$", 1.0, -0.75:0.25:0.75),
                 (:τ, "\$\\log_{10}\\,\\tau_{II}/\\tau_{bg}\$",                           1.0, -0.75:0.25:0.75),
                 (:ψ, "\$\\log_{10}\\,\\psi/\\psi_{bg}\$",                                 0.25, -0.2:0.05:0.2))
    panels  =   Any[]
    for (k, lab, lim, lev) in rows
        for (j, t) in enumerate(ratios)
            Z   =   log10.(getfield(M[t], k))
            p   =   heatmap(logΔη, θs, Z; color = cmap, clims = (-lim, lim),
                            title = k == :ε ? "a/b = $(t)" : "", xlabel = k == :ψ ? xl : "",
                            ylabel = j == 1 ? "θ [°]" : "", common...)
            contour!(p, logΔη, θs, Z; levels = collect(lev), color = :black, lw = 0.5, colorbar_entry = false)
            push!(panels, p)
        end
        tk  =   collect(range(-lim, lim, length = 5))
        push!(panels, cbar(lim, tk, tlab(tk), lab))
    end
    p1  =   plot(panels...; layout = maplayout(3, length(ratios)), size = (330 * length(ratios) + 120, 900),
                 left_margin = 8Plots.mm, bottom_margin = 3Plots.mm,
                 plot_title = "$(dname): inclusion values (colors saturate at ±1 decade; ψ at ±0.25)",
                 plot_titlefontsize = 11)
    savefig(p1, joinpath(outdir, "Maps_$(dname).png"))

    # ------------------------------------------------------------------- #
    # Figure 2: shape factor = value / value of a circle at the same Δη
    # ------------------------------------------------------------------- #
    panels  =   Any[]
    S_tk    =   [0.4, 0.5, 0.7, 1.0, 1.4, 2.0, 2.5]
    S_lev   =   [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25, 1.5, 2.0, 2.5, 3.0]
    for (j, t) in enumerate(ratios)
        Z   =   log10.(M[t].ε ./ circle.ε)
        p   =   heatmap(logΔη, θs, Z; color = cmap, clims = (-0.5, 0.5),
                        title = "a/b = $(t)", xlabel = xl, ylabel = j == 1 ? "θ [°]" : "", common...)
        contour!(p, logΔη, θs, Z; levels = log10.(S_lev), color = :black, lw = 0.5, colorbar_entry = false)
        contour!(p, logΔη, θs, Z; levels = [0.0], color = :black, lw = 1.5, colorbar_entry = false)
        push!(panels, p)
    end
    push!(panels, cbar(0.5, log10.(S_tk), tlab(S_tk), "shape factor S"))
    p2  =   plot(panels...; layout = maplayout(1, length(ratios)), size = (330 * length(ratios) + 120, 380),
                 left_margin = 8Plots.mm, bottom_margin = 6Plots.mm, top_margin = 3Plots.mm,
                 plot_title = "$(dname): shape factor S = value relative to a circle with the same Δη (same for strain rate and stress, squared for ψ; bold: S = 1)",
                 plot_titlefontsize = 10)
    savefig(p2, joinpath(outdir, "ShapeFactor_$(dname).png"))

    # ------------------------------------------------------------------- #
    # Figure 4: amplification factors A = value / background value
    # (colors and contours: A on a log scale, A = 1 (background) in bold)
    # ------------------------------------------------------------------- #
    ampt    =   vcat(1.0, ratios)
    getmap(t) = t == 1.0 ? circle : M[t]
    arows   =   (
        (k = :ε, lab = "\$A_{\\dot\\varepsilon} = \\dot{\\varepsilon}_{II}/\\dot{\\varepsilon}_{bg}\$",
         ticks = [0.2, 0.5, 1.0, 2.0, 5.0], lev = [0.25, 0.5, 0.75, 1.5, 2.0, 3.0, 4.0],
         xr = extrema(logΔη)),
        (k = :τ, lab = "\$A_{\\tau} = \\tau_{II}/\\tau_{bg}\$",
         ticks = [0.2, 0.5, 1.0, 2.0, 5.0], lev = [0.25, 0.5, 0.75, 1.5, 2.0, 3.0, 4.0],
         xr = extrema(logΔη)),
        (k = :ψ, lab = "\$A_{\\psi} = \\psi/\\psi_{bg}\$",
         ticks = [0.63, 0.8, 1.0, 1.25, 1.58], lev = [0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5],
         xr = (-2.0, 2.0)),
    )
    panels  =   Any[]
    for (r, ar) in enumerate(arows)
        lim =   log10(maximum(ar.ticks))
        for (j, t) in enumerate(ampt)
            Z   =   log10.(getfield(getmap(t), ar.k))
            p   =   heatmap(logΔη, θs, Z; color = cmap, clims = (-lim, lim), colorbar = false,
                            title = r == 1 ? (t == 1.0 ? "circle" : "a/b = $(t)") : "",
                            xlabel = r == length(arows) ? xl : "", ylabel = j == 1 ? "θ [°]" : "",
                            framestyle = :box, xlims = ar.xr, ylims = extrema(θs), yticks = 0:15:90,
                            titlefontsize = 10, guidefontsize = 9, tickfontsize = 7)
            contour!(p, logΔη, θs, Z; levels = log10.(ar.lev), color = :black, lw = 0.5, colorbar_entry = false)
            contour!(p, logΔη, θs, Z; levels = [0.0], color = :black, lw = 2.0, colorbar_entry = false)
            push!(panels, p)
        end
        push!(panels, cbar(lim, log10.(ar.ticks), tlab(ar.ticks), ar.lab))
    end
    p4  =   plot(panels...; layout = maplayout(length(arows), length(ampt)),
                 size = (300 * length(ampt) + 160, 900), left_margin = 6Plots.mm, bottom_margin = 3Plots.mm,
                 plot_title = "$(dname): amplification factors (thin contours at the levels in the text; bold: A = 1, background; ψ row zoomed to Δη = 0.01–100)",
                 plot_titlefontsize = 10)
    savefig(p4, joinpath(outdir, "Amplification_$(dname).png"))

    # ------------------------------------------------------------------- #
    # Figure 5: decomposition A_ε = C(Δη) ⋅ S(Δη, t, β) (closed form)
    # β = angle between long axis and background shortening direction
    # (pure shear: β = θ; simple shear: β = θ + 45°)
    # ------------------------------------------------------------------- #
    βs      =   collect(0.0:0.5:90.0)
    t_demo  =   5.0                                      # aspect ratio for panel (d)
    β_demo  =   [0.0, 22.5, 45.0]
    allt5   =   vcat(1.0, ratios)
    colors5 =   palette(:tab10)
    dc      =   (framestyle = :box, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8)
    qa  =   plot(10.0 .^ logΔη, IT.contrast_factor.(10.0 .^ logΔη); xscale = :log10, lw = 2,
                 color = :black, label = "", xlabel = "Δη", ylabel = "C = 2/(1+Δη)",
                 title = "(a) contrast part (circle)", ylims = (0, 2.1), dc...)
    hline!(qa, [1.0]; color = :gray, ls = :dash, label = "")
    qb  =   plot(; xlabel = "β [°]", ylabel = "\$S_0\$", xticks = 0:15:90, xlims = (0, 90),
                 title = "(b) shape part, weak limit (Δη → 0)", legend = :top, dc...)
    qc  =   plot(; xlabel = "β [°]", ylabel = "\$S_\\infty\$", xticks = 0:15:90, xlims = (0, 90),
                 title = "(c) shape part, rigid limit (Δη → ∞)", legend = false, dc...)
    for (m, t) in enumerate(allt5)
        plot!(qb, βs, [IT.shapefactor(1e-12, t, β) for β in βs]; lw = 2, color = colors5[m],
              label = t == 1.0 ? "circle" : "a/b = $(t)")
        plot!(qc, βs, [IT.shapefactor(1e12, t, β) for β in βs]; lw = 2, color = colors5[m])
    end
    qd  =   plot(; xscale = :log10, yscale = :log10, xlabel = "Δη", ylabel = "\$A_{\\dot\\varepsilon} = C \\cdot S\$",
                 title = "(d) amplification, a/b = $(t_demo)", legend = :bottomleft, dc...)
    plot!(qd, 10.0 .^ logΔη, IT.contrast_factor.(10.0 .^ logΔη); color = :black, lw = 1.5, ls = :dash,
          label = "circle (C)")
    for (m, β) in enumerate(β_demo)
        plot!(qd, 10.0 .^ logΔη, [IT.amplification(10.0^l, t_demo, β).ε for l in logΔη];
              color = colors5[m + 1], lw = 2, label = "β = $(β)°")
    end
    p5  =   plot(qa, qb, qc, qd; layout = (1, 4), size = (1600, 400), left_margin = 6Plots.mm,
                 bottom_margin = 7Plots.mm,
                 plot_title = "Amplification = contrast part C(Δη) × shape part S(Δη, a/b, β);  β: angle between long axis and shortening direction",
                 plot_titlefontsize = 11)
    savefig(p5, joinpath(outdir, "Decomposition.png"))

    # ------------------------------------------------------------------- #
    # Figure 6: A_ε and S over (β, a/b) for fixed contrasts
    # ------------------------------------------------------------------- #
    Δη_sel  =   [1e-3, 0.1, 0.5, 2.0, 10.0, 1e3]       # one column per contrast
    βg      =   collect(0.0:0.5:90.0)
    tg      =   collect(1.02:0.02:10.0)                # a/b axis (starts > 1: S ≡ 1 for a circle)
    A_tk    =   [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    A_lev   =   [0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 1.1, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]
    S6_tk   =   [0.4, 0.5, 0.7, 1.0, 1.4, 2.0, 2.5]
    S6_lev  =   [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25, 1.5, 2.0, 2.5, 3.0]
    gc      =   (framestyle = :box, xlims = (0, 90), xticks = 0:15:90, ylims = extrema(tg), colorbar = false,
                 titlefontsize = 10, guidefontsize = 9, tickfontsize = 7)
    rowA, rowS = Any[], Any[]
    for (j, Δη) in enumerate(Δη_sel)
        A   =   [IT.amplification(Δη, t, β).ε for t in tg, β in βg]
        C   =   IT.contrast_factor(Δη)
        S   =   A ./ C
        pa  =   heatmap(βg, tg, log10.(A); color = cmap, clims = (-0.8, 0.8),
                        title = @sprintf("Δη = %g,  C = %.3g", Δη, C),
                        ylabel = j == 1 ? "a/b" : "", gc...)
        contour!(pa, βg, tg, log10.(A); levels = log10.(A_lev), color = :black, lw = 0.5, colorbar_entry = false)
        minimum(A) < 1 < maximum(A) &&
            contour!(pa, βg, tg, log10.(A); levels = [0.0], color = :black, lw = 2.0, colorbar_entry = false)
        ps  =   heatmap(βg, tg, log10.(S); color = cmap, clims = (-0.5, 0.5),
                        xlabel = "β [°]", ylabel = j == 1 ? "a/b" : "", gc...)
        contour!(ps, βg, tg, log10.(S); levels = log10.(S6_lev), color = :black, lw = 0.5, colorbar_entry = false)
        contour!(ps, βg, tg, log10.(S); levels = [0.0], color = :black, lw = 2.0, colorbar_entry = false)
        push!(rowA, pa); push!(rowS, ps)
    end
    push!(rowA, cbar(0.8, log10.(A_tk), tlab(A_tk), "\$A_{\\dot\\varepsilon} = C \\cdot S\$"))
    push!(rowS, cbar(0.5, log10.(S6_tk), tlab(S6_tk), "shape factor S"))
    p6  =   plot(rowA..., rowS...; layout = maplayout(2, length(Δη_sel); wcb = 0.03), size = (290 * length(Δη_sel) + 200, 720),
                 left_margin = 8Plots.mm, right_margin = 10Plots.mm, bottom_margin = 6Plots.mm, top_margin = 2Plots.mm,
                 plot_title = "Top: amplification A_ε = C ⋅ S (bold: A = 1).  Bottom: shape factor S (bold: S = 1, neutral).  β: angle to shortening direction",
                 plot_titlefontsize = 10)
    savefig(p6, joinpath(outdir, "AmpShape_ab_beta.png"))

    # ------------------------------------------------------------------- #
    # Figure 7: A_ε over Δη for all a/b and orientations in one panel.
    # A_ε(β)² = A_ε(0°)² cos²(2β) + A_ε(45°)² sin²(2β), so for given Δη and a/b
    # every orientation lies between the β = 0° and β = 45° curves (band).
    # Panel (b): where a given β falls inside the band, S = A_ε/C vs β.
    # ------------------------------------------------------------------- #
    μb      =   10.0 .^ collect(-5.0:0.02:5.0)
    bcol    =   palette(:tab10)
    pa7 =   plot(; xscale = :log10, yscale = :log10, framestyle = :box, xlabel = "\$\\Delta\\eta = \\eta_c/\\eta_m\$",
                 ylabel = "\$A_{\\dot\\varepsilon} = \\dot{\\varepsilon}_{II}/\\dot{\\varepsilon}_{bg}\$",
                 title = "(a) band between β = 0° (dashed) and 45° (solid) = all orientations",
                 legend = :bottomleft, legendfontsize = 7, titlefontsize = 10, guidefontsize = 10, tickfontsize = 8)
    hline!(pa7, [1.0]; color = :gray, lw = 0.6, label = "")
    for (m, t) in enumerate(ratios)
        a0  =   [IT.amplification(μ, t, 0.0).ε  for μ in μb]
        a45 =   [IT.amplification(μ, t, 45.0).ε for μ in μb]
        plot!(pa7, μb, min.(a0, a45); fillrange = max.(a0, a45), fillalpha = 0.18, lw = 0,
              color = bcol[m + 1], label = "")
        plot!(pa7, μb, a45; color = bcol[m + 1], lw = 1.8, label = "a/b = $(t)")
        plot!(pa7, μb, a0;  color = bcol[m + 1], lw = 1.8, ls = :dash, label = "")
    end
    plot!(pa7, μb, IT.contrast_factor.(μb); color = :black, lw = 2, label = "circle")
    βb  =   collect(0.0:0.25:90.0)
    pb7 =   plot(; framestyle = :box, xlims = (0, 90), xticks = 0:15:90, xlabel = "β [°]",
                 ylabel = "\$S = A_{\\dot\\varepsilon}/C\$", legend = :top, legendfontsize = 7,
                 title = "(b) position inside the band: S vs β (solid: weak, dotted: rigid limit)",
                 titlefontsize = 10, guidefontsize = 10, tickfontsize = 8)
    hline!(pb7, [1.0]; color = :black, lw = 1, label = "")
    for (m, t) in enumerate(ratios)
        plot!(pb7, βb, [IT.shapefactor(1e-12, t, β) for β in βb]; color = bcol[m + 1], lw = 1.8, label = "a/b = $(t)")
        plot!(pb7, βb, [IT.shapefactor(1e12, t, β) for β in βb]; color = bcol[m + 1], lw = 1.8, ls = :dot, label = "")
    end
    p7  =   plot(pa7, pb7; layout = grid(1, 2, widths = [0.57, 0.43]), size = (1400, 540),
                 left_margin = 7Plots.mm, bottom_margin = 7Plots.mm, top_margin = 2Plots.mm,
                 plot_title = "Strain-rate amplification of an isolated elliptical inclusion (β: angle between long axis and shortening direction)",
                 plot_titlefontsize = 11)
    savefig(p7, joinpath(outdir, "AmpBand.png"))

    # ------------------------------------------------------------------- #
    # Figure 3: summary curves vs θ
    # ------------------------------------------------------------------- #
    allt    =   vcat(1.0, ratios)
    weak    =   Dict(t => [incvals(t, θ, 1e-12)[1] for θ in θs] for t in allt)
    rigid   =   Dict(t => [incvals(t, θ, 1e12)[2] for θ in θs] for t in allt)
    ψmax    =   Dict{Float64, Vector{Float64}}()
    Δηmax   =   Dict{Float64, Vector{Float64}}()
    for t in allt
        ψm, Δm  =   zeros(length(θs)), zeros(length(θs))
        for (i, θ) in enumerate(θs)
            ψv      =   [incvals(t, θ, 10.0^l)[3] for l in logΔη_fine]
            k       =   argmax(ψv)
            ψm[i], Δm[i] = ψv[k], 10.0^logΔη_fine[k]
        end
        ψmax[t], Δηmax[t] = ψm, Δm
    end
    colors  =   palette(:tab10)
    sc      =   (framestyle = :box, xlims = extrema(θs), xticks = 0:15:90, xlabel = "θ [°]",
                 titlefontsize = 10, guidefontsize = 9, tickfontsize = 8)
    pa  =   plot(; title = "weak limit (Δη → 0)", ylabel = "\$\\dot{\\varepsilon}_{II}/\\dot{\\varepsilon}_{bg}\$",
                 legend = :bottom, sc...)
    pb  =   plot(; title = "rigid limit (Δη → ∞)", ylabel = "\$\\tau_{II}/\\tau_{bg}\$", legend = false, sc...)
    pc  =   plot(; title = "maximum dissipation", ylabel = "max \$\\psi/\\psi_{bg}\$", legend = false, sc...)
    pd  =   plot(; title = "contrast at maximum dissipation", ylabel = "Δη at max ψ", yscale = :log10,
                 yticks = ([0.25, 0.5, 1.0, 2.0, 4.0], ["0.25", "0.5", "1", "2", "4"]), legend = false, sc...)
    for (m, t) in enumerate(allt)
        lab =   t == 1.0 ? "circle" : "a/b = $(t)"
        plot!(pa, θs, weak[t];  color = colors[m], lw = 2, label = lab)
        plot!(pb, θs, rigid[t]; color = colors[m], lw = 2)
        plot!(pc, θs, ψmax[t];  color = colors[m], lw = 2)
        plot!(pd, θs, Δηmax[t]; color = colors[m], lw = 2)
    end
    p3  =   plot(pa, pb, pc, pd; layout = (1, 4), size = (1500, 380), left_margin = 6Plots.mm,
                 bottom_margin = 6Plots.mm, plot_title = "$(dname): isolated elliptical inclusion",
                 plot_titlefontsize = 11)
    savefig(p3, joinpath(outdir, "Summary_$(dname).png"))

    # ------------------------------------------------------------------- #
    # Console summary and CSV output
    # ------------------------------------------------------------------- #
    println("\n$(dname): values at θ = 0° and 45° (circle: weak limit 2, rigid limit 2, max ψ = 1 at Δη = 1)")
    @printf("%6s | %18s | %18s | %24s\n", "a/b", "weak ε̇/ε̇bg 0°/45°", "rigid τ/τbg 0°/45°", "max ψ/ψbg (Δη) 0°/45°")
    i0, i45 =   findfirst(==(0.0), θs), findfirst(==(45.0), θs)
    for t in allt
        @printf("%6.1f | %8.3f / %7.3f | %8.3f / %7.3f | %5.3f (%5.2f) / %5.3f (%5.2f)\n", t,
                weak[t][i0], weak[t][i45], rigid[t][i0], rigid[t][i45],
                ψmax[t][i0], Δηmax[t][i0], ψmax[t][i45], Δηmax[t][i45])
    end
    if save_csv
        for t in ratios
            for (k, name) in ((:ε, "epsII"), (:τ, "tauII"), (:ψ, "psi"))
                IT.write_csv(joinpath(outdir, @sprintf("Map_%s_%s_ab%g.csv", dname, name, t)),
                             ["theta \\ log10_eta_ratio"; [@sprintf("%.3f", l) for l in logΔη]],
                             hcat(θs, getfield(M[t], k)))
            end
        end
        header  =   ["theta"; ["$(n)_ab$(t)" for t in allt for n in ("weak_eps", "rigid_tau", "psi_max", "eta_ratio_at_psi_max")]]
        data    =   hcat(θs, [hcat(weak[t], rigid[t], ψmax[t], Δηmax[t]) for t in allt]...)
        IT.write_csv(joinpath(outdir, "Summary_$(dname).csv"), header, data)
    end
    println("\nResults saved to $(outdir)")
    return (maps = M, circle = circle, weak = weak, rigid = rigid, ψmax = ψmax, Δηmax = Δηmax)
end

pm = run_maps();
