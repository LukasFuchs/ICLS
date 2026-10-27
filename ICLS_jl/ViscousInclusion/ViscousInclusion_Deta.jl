# ======================================================================= #
# ViscousInclusion_Deta.jl
#
# Instantaneous Stokes solution for an elliptical inclusion with given
# orientations θ under pure or simple shear, for a range of viscosity
# contrasts Δη = ηc/ηm and for one or several viscosity averaging schemes.
# For each scheme and orientation the mean, max, min, and std of ε̇_II, τ_II,
# and ψ = τ_II ⋅ ε̇_II in the inclusion and the matrix are stored (CSV)
# together with the analytical inclusion values and the errors.
# Julia counterpart of FDCSGm's ViscousInclusion_Deta.m.
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. ViscousInclusion/ViscousInclusion_Deta.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# (Re-)load the tools. Functions are called as IT.<name> instead of
# `using .InclusionTools`, so the script can be re-run in the same REPL
# session (e.g. VS Code) without name clashes between module versions.
include(joinpath(@__DIR__, "InclusionTools.jl"))
using Plots, Printf

function run_deta()
    IT = InclusionTools          # current version of the module
    # ===================== Model parameters ============================ #
    P   =   (
        L           =   1000.0,
        H           =   1000.0,
        nc          =   (x = 200, y = 200),
        a           =   300.0,
        b           =   100.0,
        ηm          =   1e23,
        ε̇bg         =   1e-15,
        deformation =   :SimpleShear,     # :PureShear or :SimpleShear
        bc          =   :analytical,
        statmask    =   :pure,          # statistics on :pure or :geometric cells
        τmethod     =   :centroid,      # τ_II = 2ηε̇_II (:centroid) or from components
    )
    avgms       =   [:geom] # [:arith, :harm, :geom]                  # averaging schemes to compare
    angles      =   [0.0, 22.5, 45.0, 90.0]                 # θ [deg]
    Δη          =   10.0 .^ range(-5, 5, length = 41)       # ηc/ηm
    outdir      =   joinpath(@__DIR__, "output", "Deta",
                        @sprintf("Ell_a_%g_b_%g_%s_nx%d", P.a, P.b, P.deformation, P.nc.x))
    # =================================================================== #
    mkpath(outdir)
    keys_   =   (:εII, :τII, :ψ)
    errkeys =   (:L2_v, :L2_v_bulk, :rms_P, :rms_P_bulk, :rms_P_mat, :rms_P_inc)
    header  =   ["eta_ratio";
                 ["$(k)_inc_$(s)" for k in keys_ for s in ("mean", "max", "min", "std")];
                 ["$(k)_mat_$(s)" for k in keys_ for s in ("mean", "max", "min", "std")];
                 ["$(k)_inc_analytical" for k in keys_];
                 [string(k) for k in errkeys]]
    col(k, s) = findfirst(==("$(k)_inc_$(s)"), header)
    cola(k)   = findfirst(==("$(k)_inc_analytical"), header)
    results =   Dict{Tuple{Symbol, Float64}, Matrix{Float64}}()

    nrun, ntot = 0, length(avgms) * length(angles) * length(Δη)
    for avgm in avgms, θ in angles
        data    =   zeros(length(Δη), length(header))
        for (i, ηr) in enumerate(Δη)
            nrun += 1
            @printf("[%4d/%d] avgm = %-5s θ = %5.1f°, ηc/ηm = %.3e\n", nrun, ntot, avgm, θ, ηr)
            r   =   IT.RunInclusionModel(; P..., avgm = avgm, θ = θ, ηc = ηr * P.ηm, verbose = false)
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
            for k in errkeys
                push!(row, r.err[k])
            end
            data[i, :] .=   row
        end
        results[(avgm, θ)] = data
        IT.write_csv(joinpath(outdir, @sprintf("Data_%s_%s_%g.csv", P.deformation, avgm, θ)),
                     header, data)
    end

    # Figures ---------------------------------------------------------- #
    ref     =   (εII = P.ε̇bg, τII = 2P.ηm * P.ε̇bg, ψ = 2P.ηm * P.ε̇bg^2)
    labels  =   (εII = "\$\\langle\\dot{\\varepsilon}_{II}\\rangle/\\dot{\\varepsilon}_{bg}\$",
                 τII = "\$\\langle\\tau_{II}\\rangle/\\tau_{bg}\$",
                 ψ   = "\$\\langle\\psi\\rangle/\\psi_{bg}\$")
    xlab    =   "\$\\Delta\\eta = \\eta_c/\\eta_m\$"
    colors  =   palette(:tab10)
    styles  =   Dict(:arith => :solid, :harm => :dash, :geom => :dot)
    ttl     =   "$(P.deformation), a/b = $(P.a/P.b), $(P.nc.x)×$(P.nc.y)"

    # (a) one figure per averaging scheme: all orientations + analytical (black, dashed)
    for avgm in avgms
        panels  =   Any[]
        for k in keys_
            p   =   plot(xscale = :log10, yscale = :log10, framestyle = :box, xlabel = xlab,
                         ylabel = labels[k], legend = k == :εII ? :bottomleft : false)
            for (j, θ) in enumerate(angles)
                d   =   results[(avgm, θ)]
                plot!(p, d[:, 1], d[:, col(k, :mean)] ./ ref[k]; color = colors[j], lw = 2,
                      label = "θ = $(θ)°")
                plot!(p, d[:, 1], d[:, cola(k)] ./ ref[k]; color = :black, lw = 1, ls = :dash,
                      label = j == 1 ? "analytical" : "")
            end
            push!(panels, p)
        end
        pA  =   plot(panels...; layout = (3, 1), size = (650, 1000), left_margin = 5Plots.mm,
                     plot_title = "$(ttl), avgm = $(avgm)", plot_titlefontsize = 11)
        savefig(pA, joinpath(outdir, "eps_tau_psi_Deta_$(P.deformation)_$(avgm).png"))

        # common decades for both axes (aspect_ratio = :equal breaks the log ticks in GR)
        allv    =   vcat([vcat(results[(avgm, θ)][:, col(:εII, :mean)] ./ ref.εII,
                               results[(avgm, θ)][:, col(:τII, :mean)] ./ ref.τII) for θ in angles]...)
        lo, hi  =   floor(log10(minimum(allv))), ceil(log10(maximum(allv)))
        tks     =   10.0 .^ (lo:hi)
        pT  =   plot(xscale = :log10, yscale = :log10, framestyle = :box,
                     xlims = (10^lo, 10^hi), ylims = (10^lo, 10^hi), xticks = tks, yticks = tks,
                     xlabel = labels[:εII], ylabel = labels[:τII], legend = :bottomleft,
                     title = "$(ttl), avgm = $(avgm)", titlefontsize = 10, size = (650, 620))
        # isolines of the dissipative energy ψ/ψ_bg = (τ/τ_bg)(ε̇/ε̇_bg) = 10^k
        for k in (lo - hi + 1):(hi - lo - 1)
            xx  =   10.0 .^ range(lo, hi, length = 50)
            plot!(pT, xx, 10.0^k ./ xx; color = :gray80, lw = 0.6, label = k == 0 ? "ψ/ψ_bg = 10^k" : "")
        end
        for (j, θ) in enumerate(angles)
            d   =   results[(avgm, θ)]
            plot!(pT, d[:, col(:εII, :mean)] ./ ref.εII, d[:, col(:τII, :mean)] ./ ref.τII;
                  color = colors[j], lw = 2, label = "θ = $(θ)°")
            plot!(pT, d[:, cola(:εII)] ./ ref.εII, d[:, cola(:τII)] ./ ref.τII;
                  color = :black, lw = 1, ls = :dash, label = j == 1 ? "analytical" : "")
        end
        savefig(pT, joinpath(outdir, "tau_vs_eps_Deta_$(P.deformation)_$(avgm).png"))
    end

    # (b) comparison of the averaging schemes: deviation from the analytical
    #     inclusion values, rows = quantity, columns = orientation
    panels  =   Any[]
    for k in keys_, (j, θ) in enumerate(angles)
        p   =   plot(xscale = :log10, framestyle = :box, legend = (k == :εII && j == 1) ? :topleft : false,
                     title = k == :εII ? "θ = $(θ)°" : "", titlefontsize = 10,
                     xlabel = k == :ψ ? xlab : "", guidefontsize = 9, tickfontsize = 7,
                     ylabel = j == 1 ? "dev. " * labels[k] * " [%]" : "")
        for avgm in avgms
            d   =   results[(avgm, θ)]
            plot!(p, d[:, 1], (d[:, col(k, :mean)] ./ d[:, cola(k)] .- 1) .* 100; color = :black,
                  ls = styles[avgm], lw = 1.5, label = string(avgm))
        end
        hline!(p, [0.0]; color = :gray, lw = 0.5, label = "")
        push!(panels, p)
    end
    pB  =   plot(panels...; layout = (3, length(angles)), size = (300 * length(angles), 850),
                 left_margin = 5Plots.mm, bottom_margin = 3Plots.mm,
                 plot_title = "$(ttl): deviation from analytical", plot_titlefontsize = 11)
    savefig(pB, joinpath(outdir, "Deviation_avgm_Deta_$(P.deformation).png"))

    # (c) bulk errors vs Δη for the averaging schemes (θ = first angle)
    θ0  =   angles[1]
    pE  =   plot(xscale = :log10, yscale = :log10, framestyle = :box, xlabel = xlab,
                 ylabel = "error [%]", legend = :outertopright, size = (800, 450),
                 title = "$(ttl), θ = $(θ0)°", titlefontsize = 10)
    ecol(name) = findfirst(==(name), header)
    for (j, (name, lab)) in enumerate((("L2_v_bulk", "L2 v bulk [% v]"), ("rms_P_mat", "rms P matrix [% τbg]"),
                                       ("rms_P_inc", "rms P incl. [% τbg]")))
        for avgm in avgms
            d   =   results[(avgm, θ0)]
            plot!(pE, d[:, 1], max.(d[:, ecol(name)], 1e-6); color = colors[j], ls = styles[avgm],
                  lw = 1.5, label = "$(lab), $(avgm)")
        end
    end
    savefig(pE, joinpath(outdir, "Errors_avgm_Deta_$(P.deformation).png"))

    println("Results saved to $(outdir)")
    return results
end

results = run_deta();
