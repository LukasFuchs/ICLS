# Viscous Inclusion

Instantaneous Stokes solutions for a circular or elliptical viscous inclusion in a viscous matrix under pure or simple shear, solved with [GeoModBox.jl](https://github.com/GeoSci-FFM/GeoModBox.jl) and compared with the analytical solution of Schmid & Podladchikov (2003). The scripts are the Julia counterparts of the `ViscousInclusion*.m` scripts of [FDCSGm](https://github.com/LukasFuchs/FDCSGm).

## Authorship and provenance

The Julia scripts in this directory (`InclusionTools.jl`, `ViscousInclusion.jl`, `ViscousInclusion_Deta.jl`, `ViscousInclusion_Rot.jl`, `resolution_test.jl`, `compare_averaging.jl`, `check_analytical.jl`, `ParameterMaps.jl`) and this README were written by Claude (Anthropic, model Claude Opus 5.5) in October 2026, at the request and under the direction of Lukas Fuchs (Goethe University Frankfurt), who also ran all models and reviewed the intermediate results.

- The analytical solution for the elliptical inclusion was implemented from Schmid & Podladchikov (2003), eqs. 48–64. Before the Julia port it was verified against an independent Python implementation (not part of this repository); the reference values in `check_analytical.jl` come from that implementation.
- The model setup and the solver loop in `RunInclusionModel` follow the viscous inclusion example of GeoModBox.jl (`examples/StokesEquation/2D/ViscousInclusion.jl`) closely and use its functions; GeoModBox itself is not modified.
- The structure of the experiments (single model, Δη sweep, rotation sweep, inclusion statistics) follows the FDCSGm MATLAB scripts. No MATLAB code was translated line by line.
- Claude could not run Julia in its own environment; all Julia results were produced on the user's machine and checked against the analytical solution.

## Content

| File | Purpose |
|---|---|
| `InclusionTools.jl` | Module with the analytical solution for an elliptical inclusion (arbitrary aspect ratio, orientation, viscosity contrast, and combination of pure and simple shear), model setup and solution with GeoModBox, diagnostics, and plotting. |
| `ViscousInclusion.jl` | Single model: fields of $\eta$, $\psi$, $\dot{\varepsilon}_{II}$, $\tau_{II}$, and comparison of $v_x$, $v_y$, $P$ with the analytical solution. |
| `ViscousInclusion_Deta.jl` | Sweep over the viscosity contrast $\Delta\eta = \eta_c/\eta_m$ for several orientations $\theta$ and averaging schemes (`avgms`), including a figure comparing the schemes. Stores mean/max/min/std of $\dot{\varepsilon}_{II}$, $\tau_{II}$, $\psi$ in inclusion and matrix (CSV) and plots them together with the analytical inclusion values. |
| `ViscousInclusion_Rot.jl` | Sweep over the orientation $\theta$ for a fixed viscosity contrast and several averaging schemes. |
| `ParameterMaps.jl` | Analytical maps of the inclusion values ε̇_II, τ_II, ψ over orientation θ and viscosity contrast Δη for several aspect ratios, the shape factor relative to a circular inclusion, the amplification factors A = value/background (e.g. ε̇_II/ε̇_bg) with labeled contours, and summary curves (weak and rigid limits, maximum dissipation), the decomposition amplification = contrast part × shape part (closed form), and maps of A_ε and S over orientation β and aspect ratio a/b for fixed contrasts. No numerical model is solved. |
| `resolution_test.jl` | Grid convergence (e.g. 50²–400²) of the inclusion means and of the bulk velocity and pressure errors for a weak and a strong inclusion and all averaging schemes, including the observed convergence order. |
| `compare_averaging.jl` | Sensitivity of errors and inclusion statistics to the viscosity averaging (`:arith`, `:harm`, `:geom`) and to the cells used for the statistics. |
| `check_analytical.jl` | Tests of the analytical solution (reference values, circular limit vs. `Dani_Solution_vec!`, continuity across the interface, far field, rigid-inclusion limit) and of the closed-form amplification factors. |

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

## Closed-form amplification factors

From the inclusion coefficient of eq. (55) (ψ_c = L z, with B₁, B₂ from eqs. 57–58 and ρ² = (t+1)/(t−1)) the amplification factors of an isolated elliptical inclusion follow in closed form. With $\mu = \Delta\eta = \eta_c/\eta_m$, aspect ratio $t = a/b \ge 1$, and $\beta$ the angle between the long axis and the background **shortening** direction:

$$
A_{\dot\varepsilon} = \frac{\dot\varepsilon_{II,inc}}{\dot\varepsilon_{bg}} = (t+1)^2 \sqrt{\frac{\sin^2 2\beta}{D_1^2} + \frac{\cos^2 2\beta}{D_2^2}}, \qquad D_1 = \mu(t^2+1) + 2t, \quad D_2 = (t^2+1) + 2\mu t,
$$

$$
A_\tau = \frac{\tau_{II,inc}}{\tau_{bg}} = \mu A_{\dot\varepsilon}, \qquad A_\psi = \frac{\psi_{inc}}{\psi_{bg}} = \mu A_{\dot\varepsilon}^2 .
$$

- Valid for any combination of pure and simple shear: only the strain-rate part of the background flow enters, the vorticity drops out. Pure shear: $\beta = \theta$; simple shear (top to the right): $\beta = \theta + 45°$ (`shortening_angle`). In words: the shortening angle is the direction in which the background flow squeezes the material most (0° for our pure shear, 135° ≡ −45° for simple shear); β is the angle between the long axis and that direction.
- Decomposition into a **contrast part** and a **shape/orientation part**:
  $A_{\dot\varepsilon} = C(\mu)\, S(\mu, t, \beta)$, with $C = 2/(1+\mu)$ (circular inclusion, Eshelby) and $S$ the shape factor shown in `ShapeFactor_*.png`. $S = 1$ for a circle or for $\mu = 1$. $S$ goes from its weak limit $S_0(t,\beta) = \tfrac{(t+1)^2}{2}\sqrt{\sin^2 2\beta/(2t)^2 + \cos^2 2\beta/(t^2+1)^2}$ through $S = 1$ at $\mu = 1$ to its rigid limit $S_\infty(t, \beta) = S_0(t, \beta + 45°)$; exactly, $S(1/\mu, t, \beta) = S(\mu, t, \beta + 45°)$. In between, $S$ is not confined to the range of these values: near orientations with $S_0 \approx 1$ it dips below 1 at moderate contrasts (a/b = 8: down to 0.86 at $\Delta\eta \approx 0.3$ or $\approx 4$).  
  For given $\mu$ and $t$, $A_{\dot\varepsilon}(\beta)^2 = A_{\dot\varepsilon}(0°)^2 \cos^2 2\beta + A_{\dot\varepsilon}(45°)^2 \sin^2 2\beta$, so all orientations lie between the $\beta = 0°$ and $45°$ curves (`AmpBand.png`).
- For a weak structure ($\mu \ll 1$): $A_{\dot\varepsilon} \approx 2 S_0(t, \beta)$, i.e. the contrast contributes at most a factor 2; the rest comes from elongation and orientation ($S_0$ up to $(t+1)^2/(4t)$ at $\beta = 45°$, i.e. $A_{\dot\varepsilon}$ up to $(t+1)^2/(2t)$).
- Implemented in `InclusionTools.jl` (`amplification`, `contrast_factor`, `shapefactor`, `shortening_angle`), tested against the full solution in `check_analytical.jl` (300 random cases, general shear), and shown in `Decomposition.png` (`ParameterMaps.jl`).

## References

Schmid, D. W., & Podladchikov, Y. Y. (2003). Analytical solutions for deformable elliptical inclusions in general shear. *Geophysical Journal International*, 155(1), 269–288. https://doi.org/10.1046/j.1365-246X.2003.02041.x
