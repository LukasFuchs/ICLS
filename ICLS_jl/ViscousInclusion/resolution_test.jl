# ======================================================================= #
# resolution_test.jl
#
# Grid convergence of the viscous inclusion models for a weak and a strong
# inclusion and several viscosity averaging schemes: deviation of the
# inclusion means from the analytical values and bulk errors of velocity and
# pressure (matrix and inclusion interior) as a function of the grid spacing.
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. -t 4 ViscousInclusion/resolution_test.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "InclusionTools.jl"))
using Plots, Printf

function run_restest()
    IT = InclusionTools
    # ===================== Parameters ================================== #
    P   =   (
        L           =   1000.0,
        H           =   1000.0,
        a           =   300.0,          # a == b -> circle; e.g. a = 300, b = 100 for an ellipse
        b           =   100.0,
        θ           =   0.0,
        ηm          =   1e23,
        ε̇bg         =   1e-15,
        deformation =   :PureShear,
        bc          =   :analytical,
        statmask    =   :pure,
        τmethod     =   :centroid,
    )
    nres        =   [50, 100, 200, 400]         # cells per direction
    ηratios     =   [1e-2, 1e2]                 # weak and strong inclusion
    avgms       =   [:arith, :harm, :geom]
    outdir      =   joinpath(@__DIR__, "output", "Resolution")
    # =================================================================== #
    mkpath(outdir)
    quantities  =   (
        (key = :dev_εII,    label = "|dev.| \$\\langle\\dot{\\varepsilon}_{II}\\rangle\$ [%]"),
        (key = :L2_v_bulk,  label = "L2 v bulk [% of v]"),
        (key = :rms_P_mat,  label = "rms P matrix [% of \$\\tau_{bg}\$]"),
        (key = :rms_P_inc,  label = "rms P inclusion [% of \$\\tau_{bg}\$]"),
    )
    header  =   ["nx", "dx", "dev_epsII", "dev_tauII", "dev_psi", "L2_v_bulk", "rms_P_mat", "rms_P_inc"]
    res     =   Dict{Tuple{Float64, Symbol}, Matrix{Float64}}()

    for ηr in ηratios, avgm in avgms
        data    =   zeros(length(nres), length(header))
        for (i, n) in enumerate(nres)
            @printf("ηc/ηm = %.0e, avgm = %-5s, %4d × %4d ... ", ηr, avgm, n, n)
            t0  =   time()
            r   =   IT.RunInclusionModel(; P..., nc = (x = n, y = n), ηc = ηr * P.ηm,
                                         avgm = avgm, verbose = false)
            dev(k)  =   (r.stats.inc[k].mean / r.ana_inc[k] - 1) * 100
            data[i, :] .= [n, P.L / n, dev(:εII), dev(:τII), dev(:ψ),
                           r.err.L2_v_bulk, r.err.rms_P_mat, r.err.rms_P_inc]
            @printf("%.1f s\n", time() - t0)
        end
        res[(ηr, avgm)] = data
        IT.write_csv(joinpath(outdir, @sprintf("Resolution_eta%.0e_%s.csv", ηr, avgm)), header, data)
    end

    # Table with observed convergence order between the two finest grids --- #
    order(e1, e2, h1, h2) = log(abs(e1) / abs(e2)) / log(h1 / h2)
    for ηr in ηratios
        @printf("\nηc/ηm = %.0e\n", ηr)
        @printf("%-6s %6s | %10s %10s %10s %10s\n", "avgm", "nx", "dev ε̇II", "L2 v bulk", "rms P mat", "rms P inc")
        for avgm in avgms
            d   =   res[(ηr, avgm)]
            for i in axes(d, 1)
                @printf("%-6s %6d | %+9.3f%% %9.4f%% %9.4f%% %9.4f%%\n", avgm, Int(d[i, 1]), d[i, 3], d[i, 6], d[i, 7], d[i, 8])
            end
            n1, n2  =   size(d, 1) - 1, size(d, 1)
            @printf("%-6s %6s | %10.2f %10.2f %10.2f %10.2f\n", avgm, "order",
                    order(d[n1, 3], d[n2, 3], d[n1, 2], d[n2, 2]),
                    order(d[n1, 6], d[n2, 6], d[n1, 2], d[n2, 2]),
                    order(d[n1, 7], d[n2, 7], d[n1, 2], d[n2, 2]),
                    order(d[n1, 8], d[n2, 8], d[n1, 2], d[n2, 2]))
        end
    end

    # Figure: errors vs grid spacing --------------------------------------- #
    colors  =   palette(:tab10)
    styles  =   Dict(:arith => :solid, :harm => :dash, :geom => :dot)
    cols    =   Dict(:dev_εII => 3, :L2_v_bulk => 6, :rms_P_mat => 7, :rms_P_inc => 8)
    panels  =   Any[]
    for ηr in ηratios, q in quantities
        p   =   plot(xscale = :log10, yscale = :log10, framestyle = :box, xflip = true,
                     xlabel = "Δx [m]", ylabel = q.label, guidefontsize = 9,
                     title = @sprintf("ηc/ηm = %.0e", ηr), titlefontsize = 10,
                     legend = q.key == :dev_εII ? :bottomleft : false)
        for (m, avgm) in enumerate(avgms)
            d   =   res[(ηr, avgm)]
            plot!(p, d[:, 2], max.(abs.(d[:, cols[q.key]]), 1e-8); color = colors[m],
                  ls = styles[avgm], marker = :circle, markersize = 3, lw = 1.5, label = string(avgm))
        end
        # reference slope: first order
        h   =   res[(ηr, avgms[1])][:, 2]
        e0  =   abs(res[(ηr, avgms[1])][1, cols[q.key]])
        plot!(p, h, e0 .* h ./ h[1]; color = :gray, ls = :dashdot, lw = 1,
              label = q.key == :dev_εII ? "O(Δx)" : "")
        push!(panels, p)
    end
    pR  =   plot(panels...; layout = (length(ηratios), length(quantities)),
                 size = (330 * length(quantities), 330 * length(ηratios) + 40),
                 left_margin = 5Plots.mm, bottom_margin = 4Plots.mm,
                 plot_title = "$(P.deformation), a = $(P.a), b = $(P.b), θ = $(P.θ)°",
                 plot_titlefontsize = 11)
    savefig(pR, joinpath(outdir, "Resolution_test.png"))
    println("\nResults saved to $(outdir)")
    return res
end

res = run_restest();
