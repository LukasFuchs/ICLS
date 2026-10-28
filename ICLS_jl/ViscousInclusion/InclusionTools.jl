# ======================================================================= #
# InclusionTools.jl
#
# Shared tools for the viscous inclusion experiments of the ICLS project:
#
#   1. Analytical solution for a deformable elliptical inclusion in general
#      (pure + simple) shear after Schmid & Podladchikov (2003), GJI 155,
#      269-288, eqs. (48)-(64). Includes arbitrary inclusion orientation.
#   2. Setup and instantaneous Stokes solution of the inclusion problem with
#      GeoModBox.jl (staggered finite differences, markers for the phases).
#   3. Diagnostics: second invariants of strain rate and deviatoric stress,
#      "dissipative energy" ψ = τ_II ⋅ ε̇_II, inclusion/matrix statistics,
#      and errors with respect to the analytical solution.
#   4. Plotting helpers.
#
# Usage (from a script in this directory):
#
#   include(joinpath(@__DIR__, "InclusionTools.jl"))
#   IT = InclusionTools                  # then call IT.RunInclusionModel(...), ...
#
#   (Avoid `using .InclusionTools` in scripts that are re-run in the same REPL
#    session: every include creates a new module version, and `using` both
#    versions makes the exported names ambiguous.)
#
# Conventions
#   - Background flow (far field), centred at the inclusion centre (xc, yc):
#         vx = ε̇ (x - xc) + γ̇ (y - yc),     vy = -ε̇ (y - yc)
#   - flow_rates(:PureShear,  ε̇bg) -> ε̇ = -ε̇bg, γ̇ = 0      (horizontal shortening,
#                                                          as in GeoModBox / FDCSGm)
#     flow_rates(:SimpleShear, ε̇bg) -> ε̇ = 0,    γ̇ = 2ε̇bg   (top-to-the-right)
#     In both cases the second invariant of the background strain rate is ε̇bg.
#   - Inclusion orientation θ in degrees, counter-clockwise from the x-axis
#     (same as the α keyword of GeoModBox's IniTracer2D).
#   - Second invariants: ε̇_II = sqrt(½(ε̇xx² + ε̇yy²) + ε̇xy²), τ_II likewise.
# ======================================================================= #
module InclusionTools

using GeoModBox
using GeoModBox.Tracers.TwoD
using GeoModBox.MomentumEquation.TwoD
using LinearAlgebra, Statistics, Printf
using Base.Threads
using Plots

export EllipseSolution, evaluate, isinside, background_velocity, flow_rates,
       inclusion_values, SetupGrid, AnalyticalFields, RunInclusionModel,
       PlotFields, PlotComparison, ellipse_outline, write_csv,
       amplification, contrast_factor, shapefactor, shortening_angle

# ======================================================================= #
# 1. Analytical solution (Schmid & Podladchikov, 2003)
# ======================================================================= #
"""
    EllipseSolution(; a, b, θ=0.0, ηm, ηc, ε̇=0.0, γ̇=0.0, xc=0.0, yc=0.0)

Precompute the coefficients of the analytical solution for an elliptical
inclusion (semi-axes `a ≥ b`, long axis rotated by `θ` degrees
counter-clockwise, centre `(xc, yc)`) with viscosity `ηc` embedded in an
infinite matrix with viscosity `ηm`, subjected to the far-field flow
`vx = ε̇(x-xc) + γ̇(y-yc)`, `vy = -ε̇(y-yc)`.

The solution follows Schmid & Podladchikov (2003, eqs. 48-64), using the
Joukowski transform `z = R(ζ + 1/ζ)`. A circle (`a == b`) is handled as the
limit `a/b → 1`, which reproduces the circular solution (Dani_Solution_vec!)
to ~1e-10.
"""
struct EllipseSolution
    a   :: Float64
    b   :: Float64
    θ   :: Float64      # [deg], counter-clockwise
    xc  :: Float64
    yc  :: Float64
    ηm  :: Float64
    ηc  :: Float64
    ε̇   :: Float64
    γ̇   :: Float64
    μ   :: Float64      # ηc / ηm
    R   :: Float64      # length scale of the Joukowski map
    ρ   :: Float64      # inclusion radius in the ζ-plane
    B   :: ComplexF64   # rotated boundary condition term, eq. (56)
    Cφ  :: ComplexF64   # φ_m = -i/2 γ̇ (ζ+1/ζ) + Cφ/ζ                 eq. (52)
    Dψ  :: ComplexF64   # ψ_m = -B (ζ+1/ζ) + Dψ/(ζ³-ζ)                eq. (53)
    K   :: ComplexF64   # φ_c = K z                                   eq. (54)
    L   :: ComplexF64   # ψ_c = L z                                   eq. (55)
end

function EllipseSolution(; a, b, θ=0.0, ηm, ηc, ε̇=0.0, γ̇=0.0, xc=0.0, yc=0.0)
    a >= b || throw(ArgumentError("a must be the major semi-axis (a ≥ b)"))
    t   =   max(a / b, 1.0 + 1e-10)         # aspect ratio; circle as limit
    ρ   =   sqrt((t + 1) / (t - 1))         # eq. (49)
    R   =   a / (ρ + 1 / ρ)
    μ   =   ηc / ηm
    α   =   -deg2rad(θ)                     # flow inclination w.r.t. ellipse
    B   =   (2ε̇ - im * γ̇) * cis(-2α)        # eq. (56)
    ρ4, ρ8 = ρ^4, ρ^8
    B1  =    μ * ρ4 + μ + ρ4 - 1            # eqs. (57)-(61)
    B2  =    μ * ρ4 - μ + ρ4 + 1
    B3  =    μ * ρ4 - μ - ρ4 + 1
    B4  =   -μ * ρ4 - μ - ρ4 + 1
    B5  =    μ * ρ8 - μ - ρ8 + 1
    Q   =   im * imag(B) / B1 - real(B) / B2
    Cφ  =   B3 * ρ^2 * Q
    Dψ  =   B5 * Q
    K   =   im * μ * B4 * γ̇ / (2B1) - ρ^2 * (μ - 1) * (im * μ * imag(B) / B1 - real(B) / B2)
    L   =   -2μ * ρ4 * (im * imag(B) / B1 + real(B) / B2)
    return EllipseSolution(a, b, θ, xc, yc, ηm, ηc, ε̇, γ̇, μ, R, ρ, B, Cφ, Dψ, K, L)
end

"""
    isinside(s::EllipseSolution, x, y)

`true` if the point `(x, y)` lies inside (or on) the elliptical inclusion.
"""
function isinside(s::EllipseSolution, x::Real, y::Real)
    sθ, cθ  =   sincosd(s.θ)
    xe      =    (x - s.xc) * cθ + (y - s.yc) * sθ
    ye      =   -(x - s.xc) * sθ + (y - s.yc) * cθ
    return (xe / s.a)^2 + (ye / s.b)^2 <= 1.0
end

"""
    evaluate(s::EllipseSolution, x, y) -> (vx, vy, P)

Analytical velocity and (dynamic) pressure at the point `(x, y)`.
Compressive pressure perturbations are positive.
"""
function evaluate(s::EllipseSolution, x::Real, y::Real)
    # rotate into the ellipse frame and normalize by R
    z   =   complex(x - s.xc, y - s.yc) * cis(-deg2rad(s.θ)) / s.R
    if isinside(s, x, y)
        # inclusion: homogeneous strain rate and stress (Eshelby)
        v   =   (s.K * z - z * conj(s.K) - conj(s.L * z)) / (2s.μ)
        p   =   -2real(s.K)
    else
        # matrix: inverse Joukowski map, branch with |ζ| > 1
        sq  =   sqrt(z * z - 4)
        ζ   =   (z + sq) / 2
        abs(ζ) < 1 && (ζ = (z - sq) / 2)
        iζ  =   inv(ζ)
        φ   =   -0.5im * s.γ̇ * (ζ + iζ) + s.Cφ * iζ
        dφ  =   -0.5im * s.γ̇ * (1 - iζ^2) - s.Cφ * iζ^2     # dφ/dζ
        ψ   =   -s.B * (ζ + iζ) + s.Dψ / (ζ^3 - ζ)
        dφz =   dφ / (1 - iζ^2)                             # dφ/dz
        v   =   (φ - z * conj(dφz) - conj(ψ)) / 2            # eq. (64), ηm = 1
        p   =   -2real(dφz)                                  # eq. (62)
    end
    v   *=  s.R * cis(deg2rad(s.θ))      # back to physical frame and units
    return real(v), imag(v), p * s.ηm
end

"""
    background_velocity(s::EllipseSolution, x, y) -> (vx, vy)

Unperturbed far-field velocity.
"""
background_velocity(s::EllipseSolution, x::Real, y::Real) =
    (s.ε̇ * (x - s.xc) + s.γ̇ * (y - s.yc), -s.ε̇ * (y - s.yc))

"""
    inclusion_values(s::EllipseSolution)

Homogeneous values inside the inclusion: `(εII, τII, ψ, P)` with
`ψ = τII ⋅ εII`. Pressure is the dynamic pressure (eq. 67).
"""
function inclusion_values(s::EllipseSolution)
    τII =   abs(s.L) * s.ηm             # |(σyy-σxx)/2 + iσxy|, eq. (73)
    εII =   τII / (2s.ηc)
    return (εII = εII, τII = τII, ψ = τII * εII, P = -2real(s.K) * s.ηm)
end

"""
    flow_rates(deformation, ε̇bg) -> (ε̇, γ̇)

Translate a deformation type and background strain-rate invariant into the
pure (`ε̇`) and simple (`γ̇`) shear rates of the far-field flow.
`deformation` is `:PureShear`, `:SimpleShear`, or a NamedTuple `(ε̇=..., γ̇=...)`
for general shear (then `ε̇bg` is ignored).
"""
function flow_rates(deformation, ε̇bg)
    if deformation isa NamedTuple
        return (deformation.ε̇, deformation.γ̇)
    elseif deformation == :PureShear
        return (-ε̇bg, 0.0)
    elseif deformation == :SimpleShear
        return (0.0, 2ε̇bg)
    else
        throw(ArgumentError("unknown deformation $deformation"))
    end
end

# ----------------------------------------------------------------------- #
# Closed-form amplification factors (derived from eqs. 52-61; see README)
# ----------------------------------------------------------------------- #
"""
    amplification(Δη, t, β) -> (ε, τ, ψ)

Exact amplification factors of an isolated elliptical inclusion in a linear
viscous matrix under any combination of pure and simple shear:

    A_ε = ε̇_II,inc / ε̇_bg = (t+1)² √( sin²(2β)/D₁² + cos²(2β)/D₂² ),
    D₁  = Δη (t²+1) + 2t,     D₂ = (t²+1) + 2Δη t,
    A_τ = τ_II,inc / τ_bg = Δη A_ε,     A_ψ = ψ_inc / ψ_bg = Δη A_ε²,

with Δη = ηc/ηm, aspect ratio t = a/b ≥ 1, and β [deg] the angle between the
long axis and the background shortening direction (see `shortening_angle`).
ε̇_bg is the second invariant of the background strain rate, τ_bg = 2ηm ε̇_bg,
ψ_bg = τ_bg ε̇_bg. The vorticity of the background flow does not enter.
"""
function amplification(Δη, t, β)
    t   >= 1 || throw(ArgumentError("t = a/b must be ≥ 1"))
    D1  =   Δη * (t^2 + 1) + 2t
    D2  =   (t^2 + 1) + 2Δη * t
    s2, c2 = sincosd(2β)
    Aε  =   (t + 1)^2 * sqrt(s2^2 / D1^2 + c2^2 / D2^2)
    return (ε = Aε, τ = Δη * Aε, ψ = Δη * Aε^2)
end

"""
    contrast_factor(Δη) = 2/(1+Δη)

Amplification of a circular inclusion (contrast part of A_ε).
"""
contrast_factor(Δη) = 2.0 / (1.0 + Δη)

"""
    shapefactor(Δη, t, β) = A_ε / contrast_factor(Δη)

Shape/orientation part of the amplification, i.e. the value relative to a
circle with the same Δη (identical for ε̇_II and τ_II; squared for ψ).
Bounded by the weak (Δη → 0) and rigid (Δη → ∞) limits, which are the same
function shifted by 45° in β.
"""
shapefactor(Δη, t, β) = amplification(Δη, t, β).ε / contrast_factor(Δη)

"""
    shortening_angle(ε̇, γ̇) -> φ [deg]

Direction of the background shortening axis (counter-clockwise from x) for
the far-field flow vx = ε̇ x + γ̇ y, vy = -ε̇ y. The angle used in
`amplification` is β = θ - φ (θ: orientation of the long axis).
Pure shear (ε̇ < 0): φ = 0; simple shear (γ̇ > 0): φ = -45° (≡ 135°).
"""
shortening_angle(ε̇, γ̇) = 0.5 * rad2deg(atan(γ̇, 2ε̇)) + 90.0

"""
    ellipse_outline(s::EllipseSolution; n=301) -> (x, y)

Coordinates of the inclusion boundary (for plotting).
"""
function ellipse_outline(s::EllipseSolution; n=301)
    φ       =   range(0, 2π, length=n)
    sθ, cθ  =   sincosd(s.θ)
    xe, ye  =   s.a .* cos.(φ), s.b .* sin.(φ)
    return s.xc .+ xe .* cθ .- ye .* sθ, s.yc .+ xe .* sθ .+ ye .* cθ
end

# ======================================================================= #
# 2. Grid, analytical fields on the staggered grid, boundary conditions
# ======================================================================= #
"""
    SetupGrid(M, NC) -> (NV, Δ, x, y)

Staggered-grid coordinates as used in the GeoModBox examples.
"""
function SetupGrid(M, NC)
    NV  =   (x = NC.x + 1, y = NC.y + 1)
    Δ   =   GridSpacing(x = (M.xmax - M.xmin) / NC.x, y = (M.ymax - M.ymin) / NC.y)
    x   =   (
        c   =   LinRange(M.xmin + Δ.x/2, M.xmax - Δ.x/2, NC.x),
        ce  =   LinRange(M.xmin - Δ.x/2, M.xmax + Δ.x/2, NC.x + 2),
        v   =   LinRange(M.xmin, M.xmax, NV.x),
    )
    y   =   (
        c   =   LinRange(M.ymin + Δ.y/2, M.ymax - Δ.y/2, NC.y),
        ce  =   LinRange(M.ymin - Δ.y/2, M.ymax + Δ.y/2, NC.y + 2),
        v   =   LinRange(M.ymin, M.ymax, NV.y),
    )
    x   =   merge(x, (
        c2d     =   x.c  .+ 0 * y.c',
        v2d     =   x.v  .+ 0 * y.v',
        vx2d    =   x.v  .+ 0 * y.ce',
        vy2d    =   x.ce .+ 0 * y.v',
    ))
    y   =   merge(y, (
        c2d     =   0 * x.c  .+ y.c',
        v2d     =   0 * x.v  .+ y.v',
        vx2d    =   0 * x.v  .+ y.ce',
        vy2d    =   0 * x.ce .+ y.v',
    ))
    return NV, Δ, x, y
end

"""
    AnalyticalFields(s, x, y, NC, NV; bc=:analytical)

Analytical pressure and velocities on the staggered grid plus the boundary
values in the same layout as GeoModBox's `Dani_Solution_vec!`:
`Pa` (centroids), `Vxa` (x.v × y.c), `Vya` (x.c × y.v), `Vx_W/E`, `Vy_S/N`
(normal velocities), `Vx_S/N`, `Vy_W/E` (tangential velocities at vertices).

With `bc=:background` the boundary values are taken from the unperturbed
far-field flow instead of the analytical solution (fields are unchanged).
"""
function AnalyticalFields(s::EllipseSolution, x, y, NC, NV; bc=:analytical)
    Pa  =   [evaluate(s, x.c[i], y.c[j])[3] for i = 1:NC.x, j = 1:NC.y]
    Vxa =   [evaluate(s, x.v[i], y.c[j])[1] for i = 1:NV.x, j = 1:NC.y]
    Vya =   [evaluate(s, x.c[i], y.v[j])[2] for i = 1:NC.x, j = 1:NV.y]
    if bc == :analytical
        f   =   (xx, yy) -> evaluate(s, xx, yy)[1:2]
    elseif bc == :background
        f   =   (xx, yy) -> background_velocity(s, xx, yy)
    else
        throw(ArgumentError("bc must be :analytical or :background"))
    end
    return (
        Pa      =   Pa,
        Vxa     =   Vxa,
        Vya     =   Vya,
        # normal velocities on the boundaries (conforming nodes)
        Vx_W    =   [f(x.v[1],   y.c[j])[1] for j = 1:NC.y],
        Vx_E    =   [f(x.v[end], y.c[j])[1] for j = 1:NC.y],
        Vy_S    =   [f(x.c[i], y.v[1]  )[2] for i = 1:NC.x],
        Vy_N    =   [f(x.c[i], y.v[end])[2] for i = 1:NC.x],
        # tangential velocities at the boundary vertices (non-conforming nodes)
        Vx_S    =   [f(x.v[i], y.v[1]  )[1] for i = 1:NV.x],
        Vx_N    =   [f(x.v[i], y.v[end])[1] for i = 1:NV.x],
        Vy_W    =   [f(x.v[1],   y.v[j])[2] for j = 1:NV.y],
        Vy_E    =   [f(x.v[end], y.v[j])[2] for j = 1:NV.y],
    )
end

# ======================================================================= #
# 3. Diagnostics
# ======================================================================= #
"""
    Invariants(ε, τ, ηc; τmethod=:centroid) -> (εII, τII, ψ)

Second invariants of strain rate and deviatoric stress on the centroids.
The shear strain rate (vertices) is averaged to the centroids.

`τmethod`:
  - `:centroid`   : τ_II = 2 η_c ε̇_II with the centroid viscosity (as in FDCSGm).
                    Consistent with the local viscosity; recommended for statistics.
  - `:components` : τ_II from the stress components of the solver (τ_xy uses the
                    vertex viscosity and is averaged to the centroids). Near the
                    interface this mixes in stresses of the neighbouring phase.
"""
function Invariants(ε, τ, ηc; τmethod = :centroid)
    av(A)   =   0.25 .* (A[1:end-1, 1:end-1] .+ A[2:end, 1:end-1] .+
                         A[1:end-1, 2:end]   .+ A[2:end, 2:end])
    εII     =   sqrt.(0.5 .* (ε.xx.^2 .+ ε.yy.^2) .+ av(ε.xy).^2)
    if τmethod == :centroid
        τII =   2.0 .* ηc .* εII
    elseif τmethod == :components
        τII =   sqrt.(0.5 .* (τ.xx.^2 .+ τ.yy.^2) .+ av(τ.xy).^2)
    else
        throw(ArgumentError("τmethod must be :centroid or :components"))
    end
    return εII, τII, τII .* εII
end

"""
    ellipse_radius(s::EllipseSolution, x, y)

Normalized elliptical radius `sqrt((xe/a)² + (ye/b)²)` of the point `(x, y)`
in the ellipse frame (`< 1` inside, `1` on the interface).
"""
function ellipse_radius(s::EllipseSolution, x::Real, y::Real)
    sθ, cθ  =   sincosd(s.θ)
    xe      =    (x - s.xc) * cθ + (y - s.yc) * sθ
    ye      =   -(x - s.xc) * sθ + (y - s.yc) * cθ
    return sqrt((xe / s.a)^2 + (ye / s.b)^2)
end

"""
    bulkmask(s, xg, yg, w)

`true` for grid points farther than about `w` from the inclusion boundary.
Uses `|r - 1| ⋅ b ≥ w` with the normalized elliptical radius `r`, which excludes
a band of half-width `w` along the flanks and a wider band (≈ w a/b) at the tips.
"""
bulkmask(s::EllipseSolution, xg, yg, w) =
    [abs(ellipse_radius(s, xg[i], yg[j]) - 1.0) * s.b >= w for i in eachindex(xg), j in eachindex(yg)]

fieldstats(A, mask) = (
    mean = mean(A[mask]), max = maximum(A[mask]),
    min  = minimum(A[mask]), std = std(A[mask]),
)

"""
    velerr(vx, vy, vxa, vya) -> (L2, Linf)

Velocity error of both components combined, relative to the analytical
velocity: L2 = ‖Δv‖/‖v_ana‖, L∞ = max|Δv_i| / max|v_ana,i| (in %). Unlike
per-component errors this stays well defined if one analytical component
vanishes (e.g. v_y for simple shear at Δη = 1).
"""
function velerr(vx, vy, vxa, vya)
    L2      =   sqrt(sum(abs2, vx .- vxa) + sum(abs2, vy .- vya)) /
                sqrt(sum(abs2, vxa) + sum(abs2, vya)) * 100.0
    Linf    =   max(maximum(abs.(vx .- vxa)), maximum(abs.(vy .- vya))) /
                max(maximum(abs, vxa), maximum(abs, vya)) * 100.0
    return L2, Linf
end

"""
    perr(ΔP, τbg) -> (rms, max)

Pressure error in % of the background stress τ_bg = 2 η_m ε̇_bg (well defined
even where the analytical pressure vanishes).
"""
perr(ΔP, τbg) = isempty(ΔP) ? (NaN, NaN) :
    (sqrt(mean(abs2, ΔP)) / τbg * 100.0, maximum(abs, ΔP) / τbg * 100.0)

# ======================================================================= #
# 4. Model run
# ======================================================================= #
"""
    RunInclusionModel(; kwargs...) -> NamedTuple

Instantaneous Stokes solution for an elliptical viscous inclusion at the
centre of a rectangular box, using GeoModBox.jl.

Keyword arguments (defaults in brackets):
  - `L`, `H`           : box width and height [1000.0, 1000.0] m (centred at 0)
  - `nc`               : number of cells `(x=..., y=...)` [(x=100, y=100)]
  - `a`, `b`           : semi-axes of the inclusion, `a ≥ b` [200.0, 200.0] m
  - `θ`                : orientation of the long axis in degrees, counter-clockwise [0.0]
  - `ηm`, `ηc`         : matrix and inclusion viscosity [1e23, 1e21] Pa s
  - `ε̇bg`              : background strain-rate invariant [1e-15] 1/s
  - `deformation`      : `:PureShear`, `:SimpleShear`, or `(ε̇=..., γ̇=...)` [:PureShear]
  - `bc`               : velocity boundary values, `:analytical` (exact infinite-
                         domain solution) or `:background` (unperturbed far field) [:analytical]
  - `nmx`, `nmy`       : markers per cell in x and y [5, 5]
  - `avgm`             : marker-to-grid viscosity averaging, `:arith`, `:harm`, `:geom` [:arith]
  - `statmask`         : cells used for the inclusion/matrix statistics:
                         `:pure` (only inclusion or only matrix markers, independent of
                         `avgm`) or `:geometric` (centroid inside/outside the ellipse) [:pure]
  - `τmethod`          : `:centroid` (τ_II = 2η ε̇_II) or `:components` [:centroid]
  - `band`             : half-width of the interface band (in cells) excluded from the
                         "bulk" errors [2]
  - `niter`            : maximum number of defect-correction iterations [10]
  - `rtol`             : relative residual tolerance [1e-10]
  - `verbose`          : print residuals and errors [true]
"""
function RunInclusionModel(;
        L = 1000.0, H = 1000.0, nc = (x = 100, y = 100),
        a = 200.0, b = 200.0, θ = 0.0,
        ηm = 1e23, ηc = 1e21, ε̇bg = 1e-15, deformation = :PureShear,
        bc = :analytical, nmx = 5, nmy = 5, avgm = :arith,
        statmask = :pure, τmethod = :centroid, band = 2,
        niter = 10, rtol = 1e-10, verbose = true)

    # Geometry and grid ------------------------------------------------- #
    M           =   Geometry(xmin = -L/2, xmax = L/2, ymin = -H/2, ymax = H/2)
    NC          =   (x = nc.x, y = nc.y)
    NV, Δ, x, y =   SetupGrid(M, NC)

    # Analytical solution ----------------------------------------------- #
    ε̇, γ̇        =   flow_rates(deformation, ε̇bg)
    sol         =   EllipseSolution(a = a, b = b, θ = θ, ηm = ηm, ηc = ηc,
                        ε̇ = ε̇, γ̇ = γ̇, xc = 0.0, yc = 0.0)
    Ana         =   AnalyticalFields(sol, x, y, NC, NV; bc = bc)

    # Fields ------------------------------------------------------------ #
    D   =   (
        vx      =   zeros(Float64, NV.x, NV.y + 1),
        vy      =   zeros(Float64, NV.x + 1, NV.y),
        Pt      =   zeros(Float64, NC...),
        ρ       =   zeros(Float64, NC...),         # no buoyancy (g = 0)
        ηc      =   zeros(Float64, NC...),
        η_ex    =   zeros(Float64, NC.x + 2, NC.y + 2),
        ph_ex   =   zeros(Float64, NC.x + 2, NC.y + 2),   # inclusion marker fraction
        ηv      =   zeros(Float64, NV...),
        wte     =   zeros(Float64, NC.x + 2, NC.y + 2),
        wtv     =   zeros(Float64, NV...),
    )
    ε   =   (xx = zeros(NC...), yy = zeros(NC...), xy = zeros(NV...))
    τ   =   (xx = zeros(NC...), yy = zeros(NC...), xy = zeros(NV...))
    divV    =   zeros(NC...)
    g       =   0.0

    # Velocity boundary conditions (Dirichlet on all sides) ------------- #
    VBC     =   (
        type    =   (E = :const, W = :const, S = :const, N = :const),
        val     =   (E = zeros(NV.y), W = zeros(NV.y), S = zeros(NV.x), N = zeros(NV.x),
                     vxE = zeros(NC.y), vxW = zeros(NC.y), vyS = zeros(NC.x), vyN = zeros(NC.x)),
    )
    VBC.val.S   .=  Ana.Vx_S
    VBC.val.N   .=  Ana.Vx_N
    VBC.val.vxW .=  Ana.Vx_W
    VBC.val.vxE .=  Ana.Vx_E
    VBC.val.W   .=  Ana.Vy_W
    VBC.val.E   .=  Ana.Vy_E
    VBC.val.vyS .=  Ana.Vy_S
    VBC.val.vyN .=  Ana.Vy_N

    # Markers: inclusion phase and viscosity ---------------------------- #
    phase   =   [0, 1]
    η       =   [ηm, ηc]
    nmark   =   nmx * nmy * NC.x * NC.y
    Aparam  =   :phase
    MAVG    =   (
        PC_th   =   [similar(D.wte) for _ = 1:nthreads()],
        PV_th   =   [similar(D.ηv)  for _ = 1:nthreads()],
        wte_th  =   [similar(D.wte) for _ = 1:nthreads()],
        wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],
    )
    Ma      =   IniTracer2D(Aparam, nmx, nmy, Δ, M, NC, 0, :Inclusion, phase;
                    ellA = a, ellB = b, α = θ)
    Markers2Cells(Ma, nmark, MAVG.PC_th, D.η_ex, MAVG.wte_th, D.wte, x, y, Δ, Aparam, η; avgm = avgm)
    D.ηc    .=  D.η_ex[2:end-1, 2:end-1]
    Markers2Vertices(Ma, nmark, MAVG.PV_th, D.ηv, MAVG.wtv_th, D.wtv, x, y, Δ, Aparam, η; avgm = avgm)
    # fraction of inclusion markers per cell (arithmetic, independent of avgm)
    Markers2Cells(Ma, nmark, MAVG.PC_th, D.ph_ex, MAVG.wte_th, D.wte, x, y, Δ, Aparam, Float64.(phase))
    phfrac  =   D.ph_ex[2:end-1, 2:end-1]

    # Numbering and system matrix --------------------------------------- #
    off     =   [NV.x * NC.y,
                 NV.x * NC.y + NC.x * NV.y]
    Num     =   (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y),
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y),
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y, NC...),
    )
    ndof    =   maximum(Num.Pt)
    F       =   zeros(ndof)
    δx      =   zeros(ndof)
    Fm      =   (x = zeros(NV.x, NC.y), y = zeros(NC.x, NV.y))
    FPt     =   zeros(NC...)
    K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
    Kfac    =   lu(K.cscmatrix)

    # Defect correction iterations -------------------------------------- #
    R0, RM  =   0.0, 0.0
    for iter = 1:niter
        Residuals2D!(D, VBC, ε, τ, divV, Δ, D.ηc, D.ηv, g, Fm, FPt)
        F[Num.Vx]   .=  Fm.x
        F[Num.Vy]   .=  Fm.y
        F[Num.Pt]   .=  FPt
        RM          =   norm(F) / length(F)
        iter == 1 && (R0 = max(RM, eps()))
        verbose && @printf("   MCE %2d: ||R|| = %1.4e, ||R||/||R₀|| = %1.4e\n", iter, RM, RM / R0)
        RM / R0 < rtol && break
        δx          .=  .-(Kfac \ F)
        D.vx[:, 2:end-1]    .+= δx[Num.Vx]
        D.vy[2:end-1, :]    .+= δx[Num.Vy]
        D.Pt                .+= δx[Num.Pt]
    end
    # strain rate and stress for the final velocity field
    Residuals2D!(D, VBC, ε, τ, divV, Δ, D.ηc, D.ηv, g, Fm, FPt)

    # Post-processing --------------------------------------------------- #
    vx_num  =   D.vx[:, 2:end-1]
    vy_num  =   D.vy[2:end-1, :]
    τbg     =   2ηm * ε̇bg

    # masks: bulk = farther than ±band cells from the interface
    w       =   band * max(Δ.x, Δ.y)
    bvx     =   bulkmask(sol, x.v, y.c, w)
    bvy     =   bulkmask(sol, x.c, y.v, w)
    bP      =   bulkmask(sol, x.c, y.c, w)
    geoP    =   [isinside(sol, x.c[i], y.c[j]) for i = 1:NC.x, j = 1:NC.y]
    bPm     =   bP .& .!geoP                    # matrix bulk
    bPi     =   bP .&   geoP                    # inclusion interior

    # The pressure is defined up to a constant. It is fixed in the matrix bulk,
    # where the solution is most accurate, and used for all pressure errors
    # and plots. All pressure errors are given in % of τ_bg.
    cP      =   mean(D.Pt[bPm] .- Ana.Pa[bPm])
    Pnum    =   D.Pt .- cP
    Pana    =   copy(Ana.Pa)
    ΔP      =   Pnum .- Pana

    L2v,  Lifv  =   velerr(vx_num, vy_num, Ana.Vxa, Ana.Vya)
    L2vb, Lifvb =   velerr(vx_num[bvx], vy_num[bvy], Ana.Vxa[bvx], Ana.Vya[bvy])
    rmsP,  maxP     =   perr(vec(ΔP), τbg)
    rmsPb, maxPb    =   perr(ΔP[bP],  τbg)
    rmsPm, maxPm    =   perr(ΔP[bPm], τbg)
    rmsPi, maxPi    =   perr(ΔP[bPi], τbg)

    err     =   (
        L2_v        = L2v,   Linf_v      = Lifv,      # [% of analytical velocity]
        L2_v_bulk   = L2vb,  Linf_v_bulk = Lifvb,
        rms_P       = rmsP,  max_P       = maxP,      # [% of τ_bg]
        rms_P_bulk  = rmsPb, max_P_bulk  = maxPb,
        rms_P_mat   = rmsPm, max_P_mat   = maxPm,
        rms_P_inc   = rmsPi, max_P_inc   = maxPi,
    )

    εII, τII, ψ =   Invariants(ε, τ, D.ηc; τmethod = τmethod)
    geomask =   [isinside(sol, x.c[i], y.c[j]) for i = 1:NC.x, j = 1:NC.y]
    if statmask == :pure
        tol     =   1e-10
        incmask =   phfrac .>= 1.0 - tol
        matmask =   phfrac .<= tol
        if !any(incmask)
            @warn "No pure inclusion cells (inclusion too small for the grid); using the geometric mask."
            incmask, matmask = geomask, .!geomask
        end
    elseif statmask == :geometric
        incmask, matmask = geomask, .!geomask
    else
        throw(ArgumentError("statmask must be :pure or :geometric"))
    end
    stats   =   (
        inc = (εII = fieldstats(εII, incmask), τII = fieldstats(τII, incmask), ψ = fieldstats(ψ, incmask)),
        mat = (εII = fieldstats(εII, matmask), τII = fieldstats(τII, matmask), ψ = fieldstats(ψ, matmask)),
    )
    ana_inc =   inclusion_values(sol)

    # background reference values
    ref     =   (εII = ε̇bg, τII = τbg, ψ = τbg * ε̇bg)

    if verbose
        @printf("   Velocity error [%% of v_ana]   all: L2 %8.4f, L∞ %8.4f | bulk (±%d cells): L2 %8.4f, L∞ %8.4f\n",
                err.L2_v, err.Linf_v, band, err.L2_v_bulk, err.Linf_v_bulk)
        @printf("   Pressure error [%% of τ_bg]    all: rms %8.4f, max %8.4f | bulk: rms %8.4f, max %8.4f\n",
                err.rms_P, err.max_P, err.rms_P_bulk, err.max_P_bulk)
        @printf("                                 matrix bulk: rms %8.4f, max %8.4f | inclusion: rms %8.4f, max %8.4f\n",
                err.rms_P_mat, err.max_P_mat, err.rms_P_inc, err.max_P_inc)
        @printf("   Statistics: %s mask (%d inclusion cells, %d of %d inside the ellipse), τ_II: %s\n",
                statmask, count(incmask), count(geomask .& incmask), count(geomask), τmethod)
        @printf("   Inclusion ⟨ε̇II⟩/ε̇bg = %.4e (analytical %.4e)\n",
                stats.inc.εII.mean / ref.εII, ana_inc.εII / ref.εII)
        @printf("   Inclusion ⟨τII⟩/τbg = %.4e (analytical %.4e)\n",
                stats.inc.τII.mean / ref.τII, ana_inc.τII / ref.τII)
        @printf("   Inclusion ⟨ψ⟩/ψbg   = %.4e (analytical %.4e)\n",
                stats.inc.ψ.mean / ref.ψ, ana_inc.ψ / ref.ψ)
        bc == :background && println("   Note: bc = :background -> errors include finite-box effects.")
    end

    return (
        params  =   (L = L, H = H, NC = NC, a = a, b = b, θ = θ, ηm = ηm, ηc = ηc,
                     ε̇bg = ε̇bg, deformation = deformation, bc = bc, avgm = avgm,
                     statmask = statmask, τmethod = τmethod, band = band),
        M = M, NC = NC, NV = NV, Δ = Δ, x = x, y = y, sol = sol, Ana = Ana,
        ηc = copy(D.ηc), ηv = copy(D.ηv),
        vx = vx_num, vy = vy_num, P = Pnum, Pana = Pana,
        εII = εII, τII = τII, ψ = ψ,
        incmask = incmask, matmask = matmask, geomask = geomask, phfrac = phfrac,
        bulkP = bP, stats = stats, ana_inc = ana_inc, ref = ref,
        err = err, residual = RM / R0,
    )
end

# ======================================================================= #
# 5. Output helpers
# ======================================================================= #
"""
    write_csv(filename, header::Vector{String}, data::AbstractMatrix)

Minimal CSV writer (no extra package required).
"""
function write_csv(filename, header, data)
    mkpath(dirname(filename))
    open(filename, "w") do io
        println(io, join(header, ","))
        for i in axes(data, 1)
            println(io, join((@sprintf("%.10e", v) for v in data[i, :]), ","))
        end
    end
    return filename
end

# ======================================================================= #
# 6. Plotting
# ======================================================================= #
_safelog(A) = log10.(max.(A, floatmin(Float64)))

"""
    PlotFields(r; quiver_step=nothing)

2×2 figure of log10(η) with velocity arrows, log10(ψ), log10(ε̇_II), and
log10(τ_II) (cf. Figs. 1 and 3 of the FDCSGm viscous inclusion README).
"""
function PlotFields(r; quiver_step = nothing)
    x, y, M =   r.x, r.y, r.M
    xo, yo  =   ellipse_outline(r.sol)
    # centroid velocities for arrows
    vxc     =   0.5 .* (r.vx[1:end-1, :] .+ r.vx[2:end, :])
    vyc     =   0.5 .* (r.vy[:, 1:end-1] .+ r.vy[:, 2:end])
    st      =   isnothing(quiver_step) ? max(1, r.NC.x ÷ 12) : quiver_step
    i0      =   max(1, st ÷ 2)
    ix, iy  =   i0:st:r.NC.x, i0:st:r.NC.y
    Xq      =   vec(x.c2d[ix, iy]);  Yq = vec(y.c2d[ix, iy])
    U, V    =   vec(vxc[ix, iy]),    vec(vyc[ix, iy])
    sc      =   0.6 * st * r.Δ.x / max(maximum(hypot.(U, V)), eps())

    common  =   (aspect_ratio = :equal, xlims = (M.xmin, M.xmax), ylims = (M.ymin, M.ymax),
                 framestyle = :box, xlabel = "x [m]", ylabel = "y [m]",
                 titlefontsize = 11, tickfontsize = 8)
    p1  =   heatmap(x.c, y.c, log10.(r.ηc)'; color = cgrad(:lapaz, rev = true),
                    title = "\$\\log_{10}\\,\\eta\$ [Pa s]", common...)
    # velocity "arrows": segments scaled with |v| and a dot at the tip
    # (fixed-size arrow heads dominate where the velocity is small)
    xs  =   vec(permutedims(hcat(Xq, Xq .+ sc .* U, fill(NaN, length(Xq)))))
    ys  =   vec(permutedims(hcat(Yq, Yq .+ sc .* V, fill(NaN, length(Yq)))))
    plot!(p1, xs, ys; color = :darkorange, lw = 1.0, label = "")
    scatter!(p1, Xq .+ sc .* U, Yq .+ sc .* V; color = :darkorange, markersize = 2,
             markerstrokewidth = 0, label = "")
    p2  =   heatmap(x.c, y.c, _safelog(r.ψ)'; color = :imola,
                    title = "\$\\log_{10}\\,\\psi = \\log_{10}(\\tau_{II}\\dot{\\varepsilon}_{II})\$", common...)
    p3  =   heatmap(x.c, y.c, _safelog(r.εII)'; color = :batlow,
                    title = "\$\\log_{10}\\,\\dot{\\varepsilon}_{II}\$ [1/s]", common...)
    p4  =   heatmap(x.c, y.c, _safelog(r.τII)'; color = :nuuk,
                    title = "\$\\log_{10}\\,\\tau_{II}\$ [Pa]", common...)
    for p in (p2, p3, p4)
        plot!(p, xo, yo; color = :black, lw = 1, ls = :dash, label = "")
    end
    return plot(p1, p2, p3, p4; layout = (2, 2), size = (1150, 1000),
                left_margin = 5Plots.mm, right_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
end

"""
    PlotComparison(r)

3×3 figure: numerical, analytical, and absolute error (in % of the maximum
analytical value) for vx, vy, and P. Velocities are normalized by
V_ref = ε̇_bg L, pressure by τ_bg = 2 η_m ε̇_bg.
"""
function PlotComparison(r)
    x, y, M =   r.x, r.y, r.M
    xo, yo  =   ellipse_outline(r.sol)
    Vref    =   r.ref.εII * (M.xmax - M.xmin)
    Pref    =   r.ref.τII
    # error normalization (consistent with r.err): velocities in % of the
    # maximum analytical velocity (both components), pressure in % of τ_bg
    vmax    =   max(maximum(abs, r.Ana.Vxa), maximum(abs, r.Ana.Vya)) / Vref
    rows    =   (
        (num = r.vx ./ Vref, ana = r.Ana.Vxa ./ Vref, xg = x.v, yg = y.c, enorm = vmax,
         name = "v_x/V_{ref}", ename = "v_x", eunit = "% of max |v|", cmap = :berlin),
        (num = r.vy ./ Vref, ana = r.Ana.Vya ./ Vref, xg = x.c, yg = y.v, enorm = vmax,
         name = "v_y/V_{ref}", ename = "v_y", eunit = "% of max |v|", cmap = :berlin),
        (num = r.P ./ Pref,  ana = r.Pana ./ Pref,    xg = x.c, yg = y.c, enorm = 1.0,
         name = "P/\\tau_{bg}", ename = "P", eunit = "% of \$\\tau_{bg}\$", cmap = :vik),
    )
    common  =   (aspect_ratio = :equal, xlims = (M.xmin, M.xmax), ylims = (M.ymin, M.ymax),
                 framestyle = :box, titlefontsize = 11, tickfontsize = 8)
    panels  =   Any[]
    for (k, rw) in enumerate(rows)
        lim     =   max(maximum(abs, rw.num), maximum(abs, rw.ana), eps())
        errf    =   abs.(rw.num .- rw.ana) ./ rw.enorm .* 100.0
        xl      =   k == 3 ? "x [m]" : ""
        p1  =   heatmap(rw.xg, rw.yg, rw.num'; color = rw.cmap, clims = (-lim, lim),
                        title = "\$$(rw.name)\$ (numerical)", xlabel = xl, ylabel = "y [m]", common...)
        p2  =   heatmap(rw.xg, rw.yg, rw.ana'; color = rw.cmap, clims = (-lim, lim),
                        title = "\$$(rw.name)\$ (analytical)", xlabel = xl, common...)
        p3  =   heatmap(rw.xg, rw.yg, errf'; color = :thermal, clims = (0.0, max(maximum(errf), eps())),
                        title = "\$$(rw.ename)\$ error [$(rw.eunit)]", xlabel = xl, common...)
        for p in (p1, p2, p3)
            plot!(p, xo, yo; color = :black, lw = 1, ls = :dash, label = "")
        end
        append!(panels, (p1, p2, p3))
    end
    return plot(panels...; layout = (3, 3), size = (1350, 1100),
                left_margin = 6Plots.mm, right_margin = 4Plots.mm, bottom_margin = 3Plots.mm)
end

end # module
