# Viscous Inclusion

Instantaneous Stokes solutions for a circular or elliptical viscous inclusion in a viscous matrix under pure or simple shear, solved with [GeoModBox.jl](https://github.com/GeoSci-FFM/GeoModBox.jl) and compared with the analytical solution of Schmid & Podladchikov (2003). The scripts are the Julia counterparts of the `ViscousInclusion*.m` scripts of [FDCSGm](https://github.com/LukasFuchs/FDCSGm).

## Content

| File | Purpose |
|---|---|
| `InclusionTools.jl` | Module with the analytical solution for an elliptical inclusion (arbitrary aspect ratio, orientation, viscosity contrast, and combination of pure and simple shear), model setup and solution with GeoModBox, diagnostics, and plotting. |
| `ViscousInclusion.jl` | Single model: fields of $\eta$, $\psi$, $\dot{\varepsilon}_{II}$, $\tau_{II}$, and comparison of $v_x$, $v_y$, $P$ with the analytical solution. |
| `ViscousInclusion_Deta.jl` | Sweep over the viscosity contrast $\Delta\eta = \eta_c/\eta_m$ for several orientations $\theta$ and averaging schemes (`avgms`), including a figure comparing the schemes. Stores mean/max/min/std of $\dot{\varepsilon}_{II}$, $\tau_{II}$, $\psi$ in inclusion and matrix (CSV) and plots them together with the analytical inclusion values. |
| `ViscousInclusion_Rot.jl` | Sweep over the orientation $\theta$ for a fixed viscosity contrast and several averaging schemes. |
| `ParameterMaps.jl` | Analytical maps of the inclusion values ε̇_II, τ_II, ψ over orientation θ and viscosity contrast Δη for several aspect ratios, the shape factor relative to a circular inclusion, and summary curves (weak and rigid limits, maximum dissipation). No numerical model is solved. |
| `resolution_test.jl` | Grid convergence (e.g. 50²–400²) of the inclusion means and of the bulk velocity and pressure errors for a weak and a strong inclusion and all averaging schemes, including the observed convergence order. |
| `compare_averaging.jl` | Sensitivity of errors and inclusion statistics to the viscosity averaging (`:arith`, `:harm`, `:geom`) and to the cells used for the statistics. |
| `check_analytical.jl` | Tests of the analytical solution (reference values, circular limit vs. `Dani_Solution_vec!`, continuity across the interface, far field, rigid-inclusion limit). |

## Usage

From the `ICLS/ICLS_jl` directory:

```bash
julia --project=. ViscousInclusion/check_analytical.jl     # run once to verify the analytical solution
julia --project=. ViscousInclusion/ViscousInclusion.jl
julia --project=. ViscousInclusion/ViscousInclusion_Deta.jl
julia --project=. ViscousInclusion/ViscousInclusion_Rot.jl
```

For multithreaded marker interpolation add, e.g., `-t 4`. The model parameters are defined at the top of each script. Results are written to `ViscousInclusion/output/` (not tracked by git).

## Setup and conventions

- Box of size $L \times H$ centred at the origin, inclusion at the centre. Phases are defined with markers (`IniTracer2D(..., :Inclusion, ...; ellA, ellB, α)`), and the viscosity is interpolated to centroids and vertices (`avgm = :arith`, `:harm`, or `:geom`).
- Velocity boundary conditions (Dirichlet on all sides):
  - `bc = :analytical`: values from the analytical (infinite-domain) solution, also for elliptical inclusions. This removes finite-box effects and allows a direct error estimate.
  - `bc = :background`: unperturbed far-field flow, i.e. the inclusion feels the box boundaries.
- Background flow, centred at the inclusion: $v_x = \dot{\varepsilon}(x-x_c) + \dot{\gamma}(y-y_c)$, $v_y = -\dot{\varepsilon}(y-y_c)$.
  - `:PureShear`: $\dot{\varepsilon} = -\dot{\varepsilon}_{bg}$ (horizontal shortening, as in GeoModBox and FDCSGm).
  - `:SimpleShear`: $\dot{\gamma} = 2\dot{\varepsilon}_{bg}$ (top to the right).
  - In both cases the second invariant of the background strain rate equals $\dot{\varepsilon}_{bg}$, so pure and simple shear results are directly comparable. Note that GeoModBox's `IniVelocity!(:SimpleShear, ...; ε)` uses $\dot{\gamma} = \varepsilon$ instead.
  - General shear: `deformation = (ε̇ = ..., γ̇ = ...)`.
- Orientation $\theta$ of the long axis in degrees, counter-clockwise from the $x$-axis.
- Second invariants $\dot{\varepsilon}_{II} = \sqrt{\tfrac{1}{2}(\dot{\varepsilon}_{xx}^2 + \dot{\varepsilon}_{yy}^2) + \dot{\varepsilon}_{xy}^2}$ (same for $\tau_{II}$) on the centroids, with the shear components averaged from the vertices. "Dissipative energy" $\psi = \tau_{II}\,\dot{\varepsilon}_{II}$.
- Inclusion and matrix statistics (`statmask`):
  - `:pure` (default): only cells without markers of the other phase (inclusion marker fraction = 1 or 0). Cells cut by the interface carry an averaged viscosity, and their stress would otherwise dominate the mean of a weak inclusion. This is independent of `avgm`.
  - `:geometric`: all centroids inside/outside the analytical ellipse.
- τ_II (`τmethod`): `:centroid` (default) uses $\tau_{II} = 2\eta_c\dot{\varepsilon}_{II}$ with the centroid viscosity, as in FDCSGm; `:components` uses the solver's stress components (τ_xy from the vertex viscosity, averaged to the centroids).
- Errors (`r.err`) are given for the whole domain and for the "bulk", i.e. excluding a band of ±`band` cells (default 2) around the interface, where the error of the sharp viscosity jump is concentrated:
  - velocity: both components combined, $\lVert \Delta \mathbf{v} \rVert / \lVert \mathbf{v}_{ana} \rVert$ (`L2_v`, `L2_v_bulk`, `Linf_*`), in %. Per-component errors are not used, since one analytical component can vanish (e.g. $v_y$ in simple shear at $\Delta\eta = 1$).
  - pressure: rms and maximum of $\Delta P$ in % of $\tau_{bg}$ (`rms_P`, `rms_P_bulk`, `rms_P_mat` for the matrix bulk, `rms_P_inc` for the inclusion interior, and `max_*`), since the analytical pressure perturbation can be zero. The pressure constant is fixed in the matrix bulk.
- Viscosity averaging: the averaging scheme is fixed within a sweep. The sweep scripts loop over `avgms` to quantify the sensitivity to the scheme, not to switch schemes with the viscosity contrast.
- Normalization in the figures: $\dot{\varepsilon}_{bg}$, $\tau_{bg} = 2\eta_m\dot{\varepsilon}_{bg}$, $\psi_{bg} = \tau_{bg}\dot{\varepsilon}_{bg}$.

## Analytical solution

`EllipseSolution` implements eqs. (48)–(64) of Schmid & Podladchikov (2003) using the Joukowski transform $z = R(\zeta + 1/\zeta)$. Inside the inclusion strain rate, stress, and pressure are homogeneous (`inclusion_values`), which gives exact reference values for the inclusion means. A circle is treated as the limit $a/b \to 1$ and reproduces GeoModBox's `Dani_Solution_vec!`.

```julia
include("InclusionTools.jl"); IT = InclusionTools
s = IT.EllipseSolution(a = 300.0, b = 100.0, θ = 30.0, ηm = 1e23, ηc = 1e19, ε̇ = -1e-15, γ̇ = 0.0)
vx, vy, P = IT.evaluate(s, 250.0, 40.0)
IT.inclusion_values(s)     # (εII, τII, ψ, P) inside the inclusion
```

The scripts call the module as `IT.<function>` rather than `using .InclusionTools`, so they can be re-run in the same REPL session (e.g. VS Code): each `include` creates a new version of the module, and `using` two versions makes the exported names ambiguous.

Differences to FDCSGm: there, the boundary velocities for elliptical inclusions were taken from the circular solution with radius $a$; here they come from the elliptical solution.

## References

Schmid, D. W., & Podladchikov, Y. Y. (2003). Analytical solutions for deformable elliptical inclusions in general shear. *Geophysical Journal International*, 155(1), 269–288. https://doi.org/10.1046/j.1365-246X.2003.02041.x
