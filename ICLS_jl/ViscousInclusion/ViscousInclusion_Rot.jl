# ======================================================================= #
# ViscousInclusion_Rot.jl
#
# Instantaneous Stokes solution for an elliptical inclusion with a given
# viscosity contrast under pure or simple shear for a range of orientations
# θ. Plots the inclusion means of ε̇_II, τ_II, and ψ = τ_II ⋅ ε̇_II over θ,
# together with the analytical inclusion values, and stores them as CSV.
# Julia counterpart of FDCSGm's ViscousInclusion_Rot.m.
#
# Run from the ICLS/Julia directory:
#   julia --project=. ViscousInclusion/ViscousInclusion_Rot.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "InclusionTools.jl"))
using .InclusionTools
using Plots, Printf

function main()
    # ===================== Model parameters ============================ #
    P   =   (
        L           =   1000.0,
        H           =   1000.0,
        nc          =   (x = 100, y = 100),
        a           =   200.0,
        b           =   100.0,
        ηm          =   1e23,
        ηc          =   1e19,           # inclusion viscosity [Pa s]
        ε̇bg         =   1e-15,
        deformation =   :PureShear,     # :PureShear or :SimpleShear
        bc          =   :analytical,
        avgm        =   :arith,
    )
    angles      =   collect(0.0:2.0:180.0)                  # θ [deg]
    outdir      =   joinpath(@__DIR__, "output", "Rot")
    # =================================================================== #
    mkpath(outdir)
    keys_   =   (:εII, :τII, :ψ)
    header  =   ["theta";
                 ["$(k)_inc_mean" for k in keys_];
                 ["$(k)_inc_analytical" for k in keys_];
                 ["L2_vx", "L2_vy", "L2_P"]]
    data    =   zeros(length(angles), length(header))

    for (i, θ) in enumerate(angles)
        @printf("θ = %5.1f° (%2d/%d)\n", θ, i, length(angles))
        r   =   RunInclusionModel(; P..., θ = θ, verbose = false)
        data[i, :] .= [θ;
                       [r.stats.inc[k].mean for k in keys_];
                       [r.ana_inc[k] for k in keys_];
                       [r.err.L2_vx, r.err.L2_vy, r.err.L2_P]]
    end
    tag =   @sprintf("%s_a%g_b%g_eta%.1e", P.deformation, P.a, P.b, P.ηc/P.ηm)
    write_csv(joinpath(outdir, "Data_Rot_$(tag).csv"), header, data)

    # Figure ----------------------------------------------------------- #
    ref     =   (εII = P.ε̇bg, τII = 2P.ηm * P.ε̇bg, ψ = 2P.ηm * P.ε̇bg^2)
    labels  =   (εII = "\$\\langle\\dot{\\varepsilon}_{II}\\rangle/\\dot{\\varepsilon}_{bg}\$",
                 τII = "\$\\langle\\tau_{II}\\rangle/\\tau_{bg}\$",
                 ψ   = "\$\\langle\\psi\\rangle/\\psi_{bg}\$")
    panels  =   Any[]
    for (j, k) in enumerate(keys_)
        p   =   plot(framestyle = :box, xlabel = "\$\\theta\$ [°]", ylabel = labels[k],
                     xlims = (0, 180), xticks = 0:30:180, legend = j == 1 ? :best : false)
        plot!(p, data[:, 1], data[:, 1 + j] ./ ref[k]; lw = 2, label = "numerical (mean)")
        plot!(p, data[:, 1], data[:, 4 + j] ./ ref[k]; lw = 1.5, ls = :dash, color = :black,
              label = "analytical")
        push!(panels, p)
    end
    pR  =   plot(panels...; layout = (3, 1), size = (650, 900), left_margin = 5Plots.mm,
                 plot_title = "$(P.deformation), a/b = $(P.a/P.b), ηc/ηm = $(P.ηc/P.ηm)",
                 plot_titlefontsize = 11)
    savefig(pR, joinpath(outdir, "Rotation_$(tag).png"))
    println("Results saved to $(outdir)")
    return data
end

data = main();
