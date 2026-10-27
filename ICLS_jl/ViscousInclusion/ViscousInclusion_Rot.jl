# ======================================================================= #
# ViscousInclusion_Rot.jl
#
# Instantaneous Stokes solution for an elliptical inclusion with a given
# viscosity contrast under pure or simple shear for a range of orientations
# θ and for one or several viscosity averaging schemes. Plots the inclusion
# means of ε̇_II, τ_II, and ψ = τ_II ⋅ ε̇_II over θ, together with the
# analytical inclusion values, and stores them as CSV.
# Julia counterpart of FDCSGm's ViscousInclusion_Rot.m.
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. ViscousInclusion/ViscousInclusion_Rot.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# (Re-)load the tools. Functions are called as IT.<name> instead of
# `using .InclusionTools`, so the script can be re-run in the same REPL
# session (e.g. VS Code) without name clashes between module versions.
include(joinpath(@__DIR__, "InclusionTools.jl"))
using Plots, Printf

function run_rot()
    IT = InclusionTools          # current version of the module
    # ===================== Model parameters ============================ #
    P   =   (
        L           =   1000.0,
        H           =   1000.0,
        nc          =   (x = 400, y = 400),
        a           =   200.0,
        b           =   100.0,
        ηm          =   1e23,
        ηc          =   1e19,           # inclusion viscosity [Pa s]
        ε̇bg         =   1e-15,
        deformation =   :PureShear,     # :PureShear or :SimpleShear
        bc          =   :analytical,
        statmask    =   :pure,          # statistics on :pure or :geometric cells
        τmethod     =   :centroid,      # τ_II = 2ηε̇_II (:centroid) or from components
    )
    avgms       =   [:geom] #[:arith, :harm, :geom]                  # averaging schemes to compare
    angles      =   collect(0.0:1.0:45.0)                  # θ [deg]
    outdir      =   joinpath(@__DIR__, "output", "Rot")
    # =================================================================== #
    mkpath(outdir)
    keys_   =   (:εII, :τII, :ψ)
    errkeys =   (:L2_v, :L2_v_bulk, :rms_P, :rms_P_bulk, :rms_P_mat, :rms_P_inc)
    header  =   ["theta";
                 ["$(k)_inc_mean" for k in keys_];
                 ["$(k)_inc_analytical" for k in keys_];
                 [string(k) for k in errkeys]]
    tag     =   @sprintf("%s_a%g_b%g_eta%.1e_nx%d", P.deformation, P.a, P.b, P.ηc/P.ηm, P.nc.x)
    results =   Dict{Symbol, Matrix{Float64}}()

    nrun, ntot = 0, length(avgms) * length(angles)
    for avgm in avgms
        data    =   zeros(length(angles), length(header))
        for (i, θ) in enumerate(angles)
            nrun += 1
            @printf("[%4d/%d] avgm = %-5s θ = %5.1f°\n", nrun, ntot, avgm, θ)
            r   =   IT.RunInclusionModel(; P..., avgm = avgm, θ = θ, verbose = false)
            data[i, :] .= [θ;
                           [r.stats.inc[k].mean for k in keys_];
                           [r.ana_inc[k] for k in keys_];
                           [r.err[k] for k in errkeys]]
        end
        results[avgm] = data
        IT.write_csv(joinpath(outdir, "Data_Rot_$(tag)_$(avgm).csv"), header, data)
    end

    # Figure ----------------------------------------------------------- #
    ref     =   (εII = P.ε̇bg, τII = 2P.ηm * P.ε̇bg, ψ = 2P.ηm * P.ε̇bg^2)
    labels  =   (εII = "\$\\langle\\dot{\\varepsilon}_{II}\\rangle/\\dot{\\varepsilon}_{bg}\$",
                 τII = "\$\\langle\\tau_{II}\\rangle/\\tau_{bg}\$",
                 ψ   = "\$\\langle\\psi\\rangle/\\psi_{bg}\$")
    colors  =   palette(:tab10)
    panels  =   Any[]
    for (j, k) in enumerate(keys_)
        p   =   plot(framestyle = :box, xlabel = "\$\\theta\$ [°]", ylabel = labels[k],
                     xlims = extrema(angles), xticks = 0:(maximum(angles) - minimum(angles) > 90 ? 30 : 5):180, legend = j == 1 ? :best : false)
        for (m, avgm) in enumerate(avgms)
            d   =   results[avgm]
            plot!(p, d[:, 1], d[:, 1 + j] ./ ref[k]; lw = 2, color = colors[m],
                  label = "numerical, $(avgm)")
        end
        d   =   results[avgms[1]]
        plot!(p, d[:, 1], d[:, 4 + j] ./ ref[k]; lw = 1.5, ls = :dash, color = :black,
              label = "analytical")
        push!(panels, p)
    end
    pR  =   plot(panels...; layout = (3, 1), size = (650, 900), left_margin = 5Plots.mm,
                 plot_title = "$(P.deformation), a/b = $(P.a/P.b), ηc/ηm = $(P.ηc/P.ηm), $(P.nc.x)×$(P.nc.y)",
                 plot_titlefontsize = 11)
    savefig(pR, joinpath(outdir, "Rotation_$(tag).png"))
    println("Results saved to $(outdir)")
    return results
end

results = run_rot();
