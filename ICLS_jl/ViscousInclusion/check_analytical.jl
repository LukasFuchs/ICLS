# ======================================================================= #
# check_analytical.jl
#
# Checks of the elliptical analytical solution in InclusionTools.jl:
#   1. Reference values from an independent (Python/numpy) implementation
#   2. Circular limit vs. GeoModBox's Dani_Solution_vec!
#   3. Velocity continuity across the inclusion boundary
#   4. Far-field limit
#   5. Rigid inclusion in pure shear, eq. (70) of Schmid & Podladchikov (2003)
#
# Run from the ICLS/Julia directory:
#   julia --project=. ViscousInclusion/check_analytical.jl
# ======================================================================= #
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "InclusionTools.jl"))
using .InclusionTools
using GeoModBox, Test

reldiff(a, b) = abs(a - b) / max(abs(b), 1e-300)

@testset "Elliptical inclusion analytical solution" begin

    @testset "1. Reference values (numpy implementation)" begin
        cases = [
            (EllipseSolution(a = 300.0, b = 100.0, θ = 30.0, ηm = 1e23, ηc = 1e19, ε̇ = 1e-15, γ̇ = 0.0),
             [(10.0, 5.0, 1.59158170685071e-14, -5.071580563363256e-15, 79987200.76795395),
              (250.0, 40.0, 4.638633826118045e-13, -3.5004136826017727e-15, 233176306.6801732),
              (-120.0, -300.0, -8.947007143367175e-14, 4.614480721954645e-13, -132253095.62373269),
              (400.0, 400.0, 3.7949532298195226e-13, -4.28142270002223e-13, -16426612.92821568)],
             (48873.220506992184, 2.443661025349609e-15, 79987200.76795395)),
            (EllipseSolution(a = 300.0, b = 100.0, θ = -45.0, ηm = 1e23, ηc = 1e27, ε̇ = 0.0, γ̇ = 2e-15,
                             xc = 500.0, yc = -500.0),
             [(510.0, -495.0, 5.001333111148142e-15, -9.997333777703718e-15, 266595567.4054324),
              (750.0, -460.0, 9.253773708177017e-15, -9.818309396390062e-14, -122276781.3540778),
              (380.0, -800.0, -5.143117570595657e-13, 9.86986132918844e-14, -91037196.57943708),
              (900.0, -100.0, 7.1681005697834e-13, -8.318994302166024e-14, -42217864.37694415)],
             (533244459.2567905, 2.666222296283953e-19, 266595567.4054324)),
            (EllipseSolution(a = 200.0, b = 50.0, θ = 110.0, ηm = 1.0, ηc = 1e-3, ε̇ = -1.0, γ̇ = 0.6),
             [(10.0, 5.0, -22.885871751687635, 40.71181101098429, 1.010056947412047),
              (250.0, 40.0, -273.6910258424224, 49.12993418267886, -0.3428220953894768),
              (-120.0, -300.0, -89.67816311944908, -367.5087660746229, 0.48261760597113207),
              (400.0, 400.0, -144.34826565842786, 417.9874514315417, 0.08194228781809776)],
             (0.005697109132823749, 2.8485545664118743, 1.010056947412047)),
        ]
        for (s, pts, inc) in cases
            for (x, y, vx, vy, p) in pts
                vxj, vyj, pj = evaluate(s, x, y)
                @test reldiff(vxj, vx) < 1e-9
                @test reldiff(vyj, vy) < 1e-9
                @test reldiff(pj, p) < 1e-9
            end
            iv = inclusion_values(s)
            @test reldiff(iv.τII, inc[1]) < 1e-9
            @test reldiff(iv.εII, inc[2]) < 1e-9
            @test reldiff(iv.P, inc[3]) < 1e-9
        end
    end

    @testset "2. Circular limit vs Dani_Solution_vec!" begin
        M       =   Geometry(xmin = 0.0, xmax = 1.0, ymin = -1.0, ymax = 0.0)
        NC      =   (x = 40, y = 40)
        NV, Δ, x, y = SetupGrid(M, NC)
        rad     =   0.2
        for (type, ε̇, γ̇) in ((:PureShear, -1.0, 0.0), (:SimpleShear, 0.0, 1.0)), ηr in (1e-4, 1e4)
            D   =   (
                Pa   = zeros(NC...), Vxa = zeros(NV.x, NC.y), Vya = zeros(NC.x, NV.y),
                Vx_N = zeros(NV.x, 1), Vx_S = zeros(NV.x, 1), Vx_W = zeros(NC.y, 1), Vx_E = zeros(NC.y, 1),
                Vy_N = zeros(NC.x, 1), Vy_S = zeros(NC.x, 1), Vy_W = zeros(NV.y, 1), Vy_E = zeros(NV.y, 1),
            )
            Dani_Solution_vec!(type, D, M, x, y, rad, ηr, NC, NV)
            s   =   EllipseSolution(a = rad, b = rad, ηm = 1.0, ηc = ηr, ε̇ = ε̇, γ̇ = γ̇,
                                    xc = 0.5, yc = -0.5)
            A   =   AnalyticalFields(s, x, y, NC, NV)
            vscale = maximum(abs, D.Vxa)
            @test maximum(abs.(A.Vxa .- D.Vxa)) / vscale < 1e-8
            @test maximum(abs.(A.Vya .- D.Vya)) / vscale < 1e-8
            @test maximum(abs.(A.Pa  .- D.Pa))  / maximum(abs, D.Pa) < 1e-8
            for f in (:Vx_W, :Vx_E, :Vy_S, :Vy_N, :Vx_S, :Vx_N, :Vy_W, :Vy_E)
                @test maximum(abs.(vec(getfield(A, f)) .- vec(getfield(D, f)))) / vscale < 1e-8
            end
        end
    end

    @testset "3. Velocity continuity across the interface" begin
        for (θ, ηr, ε̇, γ̇) in ((0.0, 1e-3, 1.0, 0.0), (30.0, 1e3, 1.0, 0.0),
                              (-50.0, 1e-2, 0.0, 1.0), (110.0, 1e4, 0.3, -1.0))
            a, b    =   0.4, 0.05
            s       =   EllipseSolution(a = a, b = b, θ = θ, ηm = 1.0, ηc = ηr, ε̇ = ε̇, γ̇ = γ̇)
            maxdiff, vmax = 0.0, 0.0
            for φ in range(0, 2π, length = 201)[1:end-1]
                xe, ye  =   a * cos(φ), b * sin(φ)
                nx, ny  =   xe / a^2, ye / b^2
                nn      =   hypot(nx, ny); nx /= nn; ny /= nn
                sθ, cθ  =   sincosd(θ)
                xb, yb  =   xe * cθ - ye * sθ, xe * sθ + ye * cθ
                nxr, nyr =  nx * cθ - ny * sθ, nx * sθ + ny * cθ
                d       =   1e-9
                vi      =   evaluate(s, xb - d * nxr, yb - d * nyr)
                vo      =   evaluate(s, xb + d * nxr, yb + d * nyr)
                maxdiff =   max(maxdiff, hypot(vi[1] - vo[1], vi[2] - vo[2]))
                vmax    =   max(vmax, hypot(vo[1], vo[2]))
            end
            @test maxdiff / vmax < 1e-6
        end
    end

    @testset "4. Far-field limit" begin
        s   =   EllipseSolution(a = 0.3, b = 0.1, θ = 20.0, ηm = 1.0, ηc = 1e-2, ε̇ = -1.0, γ̇ = 0.5)
        for (x, y) in ((200.0, -80.0), (-150.0, 170.0), (30.0, -260.0))
            vx, vy, _ = evaluate(s, x, y)
            bx, by    = background_velocity(s, x, y)
            @test hypot(vx - bx, vy - by) / hypot(bx, by) < 1e-5
        end
    end

    @testset "5. Rigid inclusion, pure shear (eq. 70)" begin
        t   =   3.0
        s   =   EllipseSolution(a = 0.3, b = 0.1, ηm = 1.0, ηc = 1e12, ε̇ = 1.0, γ̇ = 0.0)
        @test reldiff(inclusion_values(s).P / (-2.0), (t^2 - 1) / (2t)) < 1e-6
    end
end
