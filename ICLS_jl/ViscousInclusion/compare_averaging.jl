# ======================================================================= #
# compare_averaging.jl
#
# Sensitivity of errors and inclusion statistics to the marker-to-grid
# viscosity averaging (:arith, :harm, :geom) and to the choice of cells used
# for the statistics (:pure vs :geometric), for a weak and a strong inclusion.
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. ViscousInclusion/compare_averaging.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "InclusionTools.jl"))
using Printf

function run_compare()
    IT  =   InclusionTools
    P   =   (L = 1000.0, H = 1000.0, nc = (x = 200, y = 200), a = 200.0, b = 200.0, θ = 0.0,
             ηm = 1e23, ε̇bg = 1e-15, deformation = :PureShear, bc = :analytical)
    for ηr in (1e-2, 1e2)
        @printf("\nηc/ηm = %.0e  (deviation of the inclusion means from the analytical values in %%)\n", ηr)
        @printf("%-7s %-10s | %9s %9s %9s | %8s %8s %10s %10s\n", "avgm", "statmask",
                "ε̇_II", "τ_II", "ψ", "L2 v", "L2 P", "L2 v bulk", "L2 P bulk")
        for avgm in (:arith, :harm, :geom), statmask in (:pure, :geometric)
            r   =   IT.RunInclusionModel(; P..., ηc = ηr * P.ηm, avgm = avgm,
                                         statmask = statmask, verbose = false)
            dev(k)  =   (r.stats.inc[k].mean / r.ana_inc[k] - 1) * 100
            @printf("%-7s %-10s | %+8.2f%% %+8.2f%% %+8.2f%% | %7.3f%% %7.3f%% %9.3f%% %9.3f%%\n",
                    avgm, statmask, dev(:εII), dev(:τII), dev(:ψ),
                    r.err.L2_vx, r.err.L2_P, r.err.L2_vx_bulk, r.err.L2_P_bulk)
        end
    end
end

run_compare()
