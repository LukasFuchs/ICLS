# ======================================================================= #
# ViscousInclusion_Deta.jl
#
# Instantaneous Stokes solution for an elliptical inclusion with given
# orientations θ under pure or simple shear, for a range of viscosity
# contrasts Δη = ηc/ηm. For each orientation the mean, max, min, and std of
# ε̇_II, τ_II, and ψ = τ_II ⋅ ε̇_II in the inclusion and the matrix are stored
# (CSV) together with the analytical inclusion values and the errors.
# Julia counterpart of FDCSGm's ViscousInclusion_Deta.m.
#
# Run from the ICLS/Julia directory:
#   julia --project=. ViscousInclusion/ViscousInclusion_Deta.jl
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
        a           =   300.0,
        b           =   100.0,
        ηm          =   1e23,
        ε̇bg         =   1e-15,
        deformation =   :PureShear,     # :PureShear or :SimpleShear
        bc          =   :analytical,
        avgm        =   :arith,
    )
    angles      =   [0.0, 22.5, 45.0, 90.0]                 # θ [deg]
    Δη          =   10.0 .^ range(-5, 5, length = 41)       # ηc/ηm
    outdir      =   joinpath(@__DIR__, "output", "Deta",
                        @sprintf("Ell_a_%g_b_%g_%s", P.a, P.b, P.deformation))
    # =================================================================== #
    mkpath(outdir)
    keys_   =   (:εII, :τII, :ψ)
    header  =   ["eta_ratio";
                 ["$(k)_inc_$(s)" for k in keys_ for s in ("mean", "max", "min", "std")];
                 ["$(k)_mat_$(s)" for k in keys_ for s in ("mean", "max", "min", "std")];
                 ["$(k)_inc_analytical" for k in keys_];
                 ["L2_vx", "L2_vy", "L2_P"]]
    results =   Dict{Float64, Matrix{Float64}}()

    for θ in angles
        data    =   zeros(length(Δη), length(header))
        for (i, ηr) in enumerate(Δη)
            @printf("θ = %5.1f°, ηc/ηm = %.3e (%2d/%d)\n", θ, ηr, i, length(Δη))
            r   =   RunInclusionModel(; P..., θ = θ, ηc = ηr * P.ηm, verbose = false)
            row =   Float64[ηr]
            for k in keys_, s in (:mean, :max, :min, :std)
                push!(row, r.stats.inc[k][s])
            end
            for k in keys_, s in (:mean, :max, :min, :std)
                push!(row, r.stats.mat[k][s])
            end
            for k in keys_
                push!(row, r.ana_inc[k])
            end
            append!(row, (r.err.L2_vx, r.err.L2_vy, r.err.L2_P))
            data[i, :] .=   row
        end
        results[θ] = data
        write_csv(joinpath(outdir, @sprintf("Data_%s_%g.csv", P.deformation, θ)), header, data)
    end

    # Figures ---------------------------------------------------------- #
    ref     =   (εII = P.ε̇bg, τII = 2P.ηm * P.ε̇bg, ψ = 2P.ηm * P.ε̇bg^2)
    col(k, s) = findfirst(==("$(k)_inc_$(s)"), header)
    cola(k)   = findfirst(==("$(k)_inc_analytical"), header)
    labels  =   (εII = "\$\\langle\\dot{\\varepsilon}_{II}\\rangle/\\dot{\\varepsilon}_{bg}\$",
                 τII = "\$\\langle\\tau_{II}\\rangle/\\tau_{bg}\$",
                 ψ   = "\$\\langle\\psi\\rangle/\\psi_{bg}\$")
    colors  =   palette(:tab10)

    panels  =   Any[]
    for k in keys_
        p   =   plot(xscale = :log10, yscale = :log10, framestyle = :box,
                     xlabel = "\$\\Delta\\eta = \\eta_c/\\eta_m\$", ylabel = labels[k],
                     legend = k == :εII ? :bottomleft : false)
        for (j, θ) in enumerate(angles)
            d   =   results[θ]
            plot!(p, d[:, 1], d[:, col(k, :mean)] ./ ref[k]; color = colors[j], lw = 2,
                  label = "θ = $(θ)°")
            plot!(p, d[:, 1], d[:, cola(k)] ./ ref[k]; color = colors[j], lw = 1.5, ls = :dash,
                  label = "")
        end
        push!(panels, p)
    end
    p1  =   plot(panels...; layout = (3, 1), size = (650, 1000), left_margin = 5Plots.mm,
                 plot_title = "$(P.deformation), a/b = $(P.a/P.b)  (dashed: analytical)",
                 plot_titlefontsize = 11)

    p2  =   plot(xscale = :log10, yscale = :log10, framestyle = :box, aspect_ratio = :equal,
                 xlabel = labels[:εII], ylabel = labels[:τII], legend = :bottomleft,
                 title = "$(P.deformation), a/b = $(P.a/P.b)", size = (650, 600))
    for (j, θ) in enumerate(angles)
        d   =   results[θ]
        plot!(p2, d[:, col(:εII, :mean)] ./ ref.εII, d[:, col(:τII, :mean)] ./ ref.τII;
              color = colors[j], lw = 2, label = "θ = $(θ)°")
    end

    savefig(p1, joinpath(outdir, "eps_tau_psi_Deta_$(P.deformation).png"))
    savefig(p2, joinpath(outdir, "tau_vs_eps_Deta_$(P.deformation).png"))
    println("Results saved to $(outdir)")
    return results
end

results = main();
