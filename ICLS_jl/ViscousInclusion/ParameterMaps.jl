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
#   + CSV files with the maps and the summary curves.
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
    # Figure 1: absolute values, diverging colors centered at background
    # ------------------------------------------------------------------- #
    xl      =   "\$\\log_{10}\\,\\Delta\\eta\$"
    common  =   (framestyle = :box, xlims = extrema(logΔη), ylims = extrema(θs),
                 yticks = 0:15:90, titlefontsize = 10, guidefontsize = 9, tickfontsize = 7)
    rows    =   ((:ε, "\$\\log_{10}\\,\\dot{\\varepsilon}_{II}/\\dot{\\varepsilon}_{bg}\$", 1.0, -1.0:0.25:1.0),
                 (:τ, "\$\\log_{10}\\,\\tau_{II}/\\tau_{bg}\$",                           1.0, -1.0:0.25:1.0),
                 (:ψ, "\$\\log_{10}\\,\\psi/\\psi_{bg}\$",                                 0.25, [-1.0, -0.5, -0.2, -0.1, 0.0, 0.05, 0.1, 0.15, 0.2]))
    panels  =   Any[]
    for (k, lab, lim, lev) in rows, (j, t) in enumerate(ratios)
        Z   =   log10.(getfield(M[t], k))
        p   =   heatmap(logΔη, θs, Z; color = cgrad(:RdBu, rev = true), clims = (-lim, lim),
                        colorbar = j == length(ratios), colorbar_title = j == length(ratios) ? lab : "",
                        title = k == :ε ? "a/b = $(t)" : "", xlabel = k == :ψ ? xl : "",
                        ylabel = j == 1 ? "θ [°]" : "", common...)
        contour!(p, logΔη, θs, Z; levels = collect(lev), color = :black, lw = 0.5, colorbar_entry = false)
        push!(panels, p)
    end
    p1  =   plot(panels...; layout = (3, length(ratios)), size = (330 * length(ratios) + 80, 900),
                 left_margin = 4Plots.mm, bottom_margin = 3Plots.mm,
                 plot_title = "$(dname): inclusion values (colors saturate at ±1 decade; ψ at ±0.25)",
                 plot_titlefontsize = 11)
    savefig(p1, joinpath(outdir, "Maps_$(dname).png"))

    # ------------------------------------------------------------------- #
    # Figure 2: shape factor = value / value of a circle at the same Δη
    # ------------------------------------------------------------------- #
    panels  =   Any[]
    lev     =   [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]
    for (j, t) in enumerate(ratios)
        Z   =   log10.(M[t].ε ./ circle.ε)
        p   =   heatmap(logΔη, θs, Z; color = cgrad(:RdBu, rev = true), clims = (-0.5, 0.5),
                        colorbar = j == length(ratios),
                        colorbar_title = j == length(ratios) ? "\$\\log_{10}\$ shape factor" : "",
                        title = "a/b = $(t)", xlabel = xl, ylabel = j == 1 ? "θ [°]" : "", common...)
        contour!(p, logΔη, θs, Z; levels = lev, color = :black, lw = 0.5, colorbar_entry = false)
        push!(panels, p)
    end
    p2  =   plot(panels...; layout = (1, length(ratios)), size = (330 * length(ratios) + 80, 380),
                 left_margin = 8Plots.mm, bottom_margin = 6Plots.mm, top_margin = 3Plots.mm,
                 plot_title = "$(dname): shape factor = value relative to a circle with the same Δη (equal for strain rate and stress, squared for ψ)",
                 plot_titlefontsize = 10)
    savefig(p2, joinpath(outdir, "ShapeFactor_$(dname).png"))

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
                 legend = false, sc...)
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
