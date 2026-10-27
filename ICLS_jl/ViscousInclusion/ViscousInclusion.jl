# ======================================================================= #
# ViscousInclusion.jl
#
# Instantaneous Stokes solution for a single circular or elliptical viscous
# inclusion under pure or simple shear, compared with the analytical solution
# of Schmid & Podladchikov (2003). Julia counterpart of FDCSGm's
# ViscousInclusion.m (which supported the analytical comparison only for
# circles).
#
# Run from the ICLS/ICLS_jl directory:
#   julia --project=. ViscousInclusion/ViscousInclusion.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# (Re-)load the tools. Functions are called as IT.<name> instead of
# `using .InclusionTools`, so the script can be re-run in the same REPL
# session (e.g. VS Code) without name clashes between module versions.
include(joinpath(@__DIR__, "InclusionTools.jl"))
using Plots, Printf

function run_single()
    IT = InclusionTools          # current version of the module
    # ===================== Model parameters ============================ #
    P   =   (
        L           =   1000.0,         # box width [m]
        H           =   1000.0,         # box height [m]
        nc          =   (x = 200, y = 200),
        a           =   200.0,          # major semi-axis [m]
        b           =   200.0,          # minor semi-axis [m]  (a == b -> circle)
        θ           =   0.0,            # orientation [deg], counter-clockwise
        ηm          =   1e23,           # matrix viscosity [Pa s]
        ηc          =   1e21,           # inclusion viscosity [Pa s]
        ε̇bg         =   1e-15,          # background strain-rate invariant [1/s]
        deformation =   :PureShear,     # :PureShear, :SimpleShear, or (ε̇=..., γ̇=...)
        bc          =   :analytical,    # :analytical or :background
        avgm        =   :arith,         # viscosity averaging: :arith, :harm, :geom
        statmask    =   :pure,          # statistics on :pure or :geometric cells
        τmethod     =   :centroid,      # τ_II = 2ηε̇_II (:centroid) or from components
    )
    save_fig    =   true
    outdir      =   joinpath(@__DIR__, "output", "single")
    # =================================================================== #

    println("Viscous inclusion: a = $(P.a), b = $(P.b), θ = $(P.θ)°, ",
            "ηc/ηm = $(P.ηc/P.ηm), $(P.deformation), bc = $(P.bc)")
    r   =   IT.RunInclusionModel(; P...)

    # Inclusion and matrix statistics --------------------------------- #
    println("\n   Inclusion / matrix statistics (normalized by background values):")
    for (name, key) in (("ε̇_II", :εII), ("τ_II", :τII), ("ψ", :ψ))
        si, sm, ref = r.stats.inc[key], r.stats.mat[key], r.ref[key]
        @printf("   %-5s inc: mean %.4e  min %.4e  max %.4e  std %.4e | mat: mean %.4e\n",
                name, si.mean/ref, si.min/ref, si.max/ref, si.std/ref, sm.mean/ref)
    end

    # Figures --------------------------------------------------------- #
    pF  =   IT.PlotFields(r)
    pC  =   IT.PlotComparison(r)
    if save_fig
        mkpath(outdir)
        tag =   @sprintf("%s_a%g_b%g_th%g_eta%.1e", P.deformation isa Symbol ? P.deformation : "GeneralShear",
                         P.a, P.b, P.θ, P.ηc/P.ηm)
        savefig(pF, joinpath(outdir, "Fields_$(tag).png"))
        savefig(pC, joinpath(outdir, "Comparison_$(tag).png"))
        println("\n   Figures saved to $(outdir)")
    else
        display(pF); display(pC)
    end
    return r
end

r = run_single();
