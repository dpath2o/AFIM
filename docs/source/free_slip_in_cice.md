# Free-slip Boundary Condition in CICE6 (AFIM / CICE\_free-slip)

*Last updated: 2026-01-15*

This note documents the **free-slip boundary condition** implementation added to CICE6 in the AFIM/CICE\_free-slip development branch, and how it interacts with the **Coastal Drag Parameterisation (CDP)** used to simulate Antarctic landfast ice.

It is written to be placed directly in `docs/source/` (ReadTheDocs/Sphinx). It complements the form-factor documentation in `lateral_drag_form_factor_creation.md`.

---

## Contents

1. [Purpose and scope](#1-purpose-and-scope)  
2. [Background: why free-slip matters](#2-background-why-free-slip-matters)  
3. [`ice_in` configuration](#3-ice_in-configuration)  
4. [Code changes and control flow](#4-code-changes-and-control-flow)  
5. [Free-slip strain-rate discretisation](#5-free-slip-strain-rate-discretisation)  
6. [Interaction with Coastal Drag Parameterisation (CDP)](#6-interaction-with-coastal-drag-parameterisation-cdp)  
7. [Deriving and interpreting **C\_s**](#7-deriving-and-interpreting-c_s)  
8. [Resolution dependence: why **C\_s** can vary by 10–1000×](#8-resolution-dependence-why-c_s-can-vary-by-101000)  
9. [Verification checklist](#9-verification-checklist)  
10. [Caveats and known limitations](#10-caveats-and-known-limitations)  
11. [References](#11-references)  

---

## 1. Purpose and scope

CICE historically behaves as **no-slip** along coastlines in many configurations, especially on a **B-grid**, because the velocity points adjacent to land are treated as inactive (or effectively pinned). This can:

- over-damp along-coast motion,
- over-stabilise coastal ice and promote unrealistically persistent “stuck” ice in channels, and
- contribute to excessive fast-ice build-up in some geometries.

The AFIM/CICE\_free-slip branch adds a selectable **free-slip** option intended to:

- maintain *no normal flow* into land, while
- avoiding artificial suppression of *tangential* motion and shear at coastal boundaries,
- without changing the EVP rheology itself.

This implementation is designed for the **C-grid** momentum discretisation (`grid_ice = "C"` or `"CD"` variants), where normal/tangential components relative to a boundary can be represented more cleanly.

---

## 2. Background: why free-slip matters

### 2.1 No-slip vs free-slip (conceptual)

Let **u** be the ice velocity and **n** the outward unit normal to a boundary.

- **No-slip**:  

  \[
  \mathbf{u}\cdot\mathbf{n} = 0 \quad \text{and} \quad \mathbf{u}\cdot\mathbf{t} = 0
  \]

  (both normal and tangential components are forced to zero at the boundary).

- **Free-slip** (idealised):  

  \[
  \mathbf{u}\cdot\mathbf{n} = 0 \quad \text{and} \quad \frac{\partial (\mathbf{u}\cdot\mathbf{t})}{\partial n} = 0
  \]

  (no penetration, but tangential flow is permitted with zero normal gradient—i.e., no imposed boundary-layer shear).

In discrete sea-ice momentum solvers, the *practical* implementation typically happens inside the **strain-rate** calculation, where derivatives require neighbour velocities. The primary numerical distinction is how “missing” neighbour values across land are handled when forming gradients.

### 2.2 B-grid vs C-grid note

On a **B-grid**, the velocity point can sit in a way that makes coastline handling more “all-or-nothing” (velocity points can be deactivated when any corner of the cell is land), which behaves similarly to no-slip. On a **C-grid**, velocity components sit on faces, which better supports enforcing no-normal-flow while retaining tangential motion.

For additional discussion of how grid staggering affects sea-ice dynamics and coastlines, see the references in Section 11.

---

## 3. `ice_in` configuration

The free-slip behaviour is controlled through the namelist entry:

- `boundary_condition = 'free_slip'` *(default remains no-slip unless specified)*

In addition, if CDP is enabled, the following are relevant:

- `coastal_drag = .true.`  
- `Cs` *(static coefficient, units m s\ :sup:`-2`)*  
- `u0` *(residual speed, units m s\ :sup:`-1`)*  
- the **form-factor** configuration (F2 file and mapping method)

### 3.1 Example `ice_in` snippet

Below is a minimal excerpt consistent with the AFIM configuration you provided (showing both free-slip and CDP):

```fortran
&dynamics_nml
    boundary_condition    = 'free_slip'
    coastal_drag          = .true.
    Cs                   = 2e-4
    u0                   = 5e-4
/

&grid_nml
    F2_file           = '/path/to/form_factor_file.nc'
    F2x_varname       = 'F2x'
    F2y_varname       = 'F2y'
    F2_map_method     = 'max'     ! 'avg' or 'max'
    F2_test           = .false.
    F2_test_val       = 0.25
/
```

**Notes**

- `boundary_condition` is parsed as a string and used to set an internal logical `noslip` flag.
- `F2_map_method` controls how *T-grid* form factors are mapped onto **velocity faces**:
  - `'avg'`: average adjacent T-cell values to the face
  - `'max'`: take the maximum of adjacent values (more aggressive drag near rough coastlines)

---

## 4. Code changes and control flow

This section summarises *where* free-slip is applied and how it threads through the solver.

### 4.1 `ice_dyn_evp.F90`

Key additions/changes (AFIM/CICE\_free-slip branch):

1. **Initialisation**  
   In `init_evp`, after `init_dyn_shared(dt_dyn)`, form factors are loaded when CDP is enabled:

   - `if (coastal_drag) call load_F2_form_factors()`

2. **Boundary-condition selection**
   A `boundary_condition` string is interpreted and a logical `noslip` is set, e.g.:

   - `noslip = (trim(boundary_condition) == 'no_slip')`

3. **C-grid strain-rate call site**
   In the C-grid EVP pathway, strain rates at the U-point are computed via either:

   - `strain_rates_U_no_slip(...)` *(legacy behaviour)*  
   - `strain_rates_U_free_slip(...)` *(new behaviour)*  

   depending on `noslip`.

4. **CDP “Ku” update**
   After `stepu_C`/`stepv_C`, the coastal stress factor is formed on E and N faces:

   - `KuE = emass * F2E * Cs`
   - `KuN = nmass * F2N * Cs`

   implemented via `coastal_drag_stress_factor(...)`.

### 4.2 `ice_grid.F90`

A new routine `load_F2_form_factors()` reads pre-computed form factors from a NetCDF file and maps them onto:

- `F2E` *(E-faces)*
- `F2N` *(N-faces)*

Key behaviour:

- NaNs are converted to zero.
- A mapping from **T-projection** to **faces** is applied (`avg` or `max`).
- Values are scattered to blocks and halo-updated.
- A fallback **test mode** (`F2_test`) can apply uniform form factors on the domain perimeter.

### 4.3 `ice_dyn_shared.F90`

This module carries most of the free-slip implementation in practice:

- Defines the `boundary_condition` configuration variable.
- Implements `strain_rates_U_free_slip(...)`.
- Adds CDP parameters/arrays (`coastal_drag`, `Cs`, `u0`, `KuE/KuN`, `Kux/Kuy`, etc.).
- Updates `stepu_C` and `stepv_C` to include a CDP contribution in the implicit momentum solve.

---

## 5. Free-slip strain-rate discretisation

### 5.1 Where free-slip is applied numerically

In this implementation, **free-slip is enforced by how neighbour velocities are supplied to finite-difference gradients** used in the U-point strain-rate calculation. The goal is to avoid forcing tangential components to zero solely because the adjacent point lies on land.

The free-slip routine uses the existing edge masks:

- `epm` — E-face activity mask  
- `npm` — N-face activity mask

These are typically 1 for active ocean/ice velocity points and 0 for land/inactive.

### 5.2 Reflection-style neighbour substitution (key idea)

When a neighbour velocity lies across land (mask = 0), the free-slip formulation replaces it with a reflected/copy value from the interior side, e.g. schematically:

- If `epm(i,j) = 1` but `epm(i,j+1) = 0` (north neighbour is land), then  
  **use `vvelE(i,j)` as the “north neighbour”** when forming \(\partial v/\partial y\).

This is implemented by blending masked neighbours, for example (illustrative form):

\[
v^{E}_{i,j+1} \leftarrow v^{E}_{i,j+1}\,epm_{i,j+1} + (epm_{i,j}-epm_{i,j+1})\,epm_{i,j}\,v^{E}_{i,j}
\]

which reduces to:

- \(v^{E}_{i,j+1}\) when `epm(i,j+1)=1`, and
- \(v^{E}_{i,j}\) when `epm(i,j+1)=0` and `epm(i,j)=1`.

Analogous substitutions are performed for the other neighbour directions/components needed by the U-point strain-rate operator.

### 5.3 Effect on the strain-rate tensor

Using interior-value substitution across land effectively enforces:

- **no normal flow** at land boundaries (through the underlying masks/velocity placement), and
- **zero normal gradient of tangential velocity** in the discrete derivative (free-slip analogue)

rather than forcing the tangential component itself to zero.

---

## 6. Interaction with Coastal Drag Parameterisation (CDP)

Free-slip and CDP address different physical/numerical issues:

- **Free-slip** relaxes artificial tangential pinning imposed by boundary handling.
- **CDP** adds an *explicit* lateral momentum sink tied to coastline roughness and ice mass.

In AFIM/CICE\_free-slip, these can be used together: free-slip provides the numerically permissive boundary, while CDP supplies the physically motivated coastal damping needed to sustain realistic landfast ice.

### 6.1 CDP stress law (from Liu et al. 2022)

The lateral (coastal) drag stress is written as (one of the Liu et al. variants):

\[
\boldsymbol{\tau}_{\ell} = m\,F_2\,C_s\,\frac{\mathbf{u}}{\lvert\mathbf{u}\rvert + u_0}
\]

where:

- \(m\) is sea-ice mass per unit area (kg m\ :sup:`-2`), typically \(m = \rho_i\,h\,a\)
- \(F_2\) is a **dimensionless** form factor representing effective coastline roughness/orientation
- \(C_s\) is a tunable coefficient with units of **acceleration** (m s\ :sup:`-2`)
- \(u_0\) is a small residual speed (m s\ :sup:`-1`) preventing a singularity as \(|u|\to 0\)

This formulation transitions between:

- linear drag at very small speeds (regularised by \(u_0\)), and
- a capped (approximately constant-magnitude) stress as \(|u|\gg u_0\).

### 6.2 How CDP enters the C-grid momentum solve

Implementation pattern:

1. Compute a stress scale on faces:
   \[
   K_u = m\,F_2\,C_s
   \]
   where `m` is taken as `emass` or `nmass` (face mass per area).

2. Convert to an implicit drag coefficient:
   \[
   C_\ell = \frac{K_u}{\lvert\mathbf{u}\rvert + u_0}
   \]
   (units kg m\ :sup:`-2` s\ :sup:`-1`).

3. Add \(C_\ell\) into the implicit coefficient `cca` in `stepu_C`/`stepv_C`.

4. Export diagnostic stress components:
   - `Kux = -u * C_ell`
   - `Kuy = -v * C_ell`

These outputs are useful for confirming where and how strongly CDP is acting.

---

## 7. Deriving and interpreting **C\_s**

### 7.1 Units and meaning

From the stress scale:

\[
K_u = m\,F_2\,C_s
\]

- \(m\) has units kg m\ :sup:`-2`
- \(F_2\) is dimensionless
- \(C_s\) has units m s\ :sup:`-2`

Therefore:

- \(K_u\) has units kg m\ :sup:`-1` s\ :sup:`-2` = N m\ :sup:`-2` = Pa

So **C\_s is an acceleration scale**. You can interpret it as:

> the maximum coastal-drag deceleration per unit ice mass, in a cell with \(F_2=1\), when the stress is fully “activated”.

A useful back-of-envelope conversion:

- \(C_s = 10^{-4}\,\mathrm{m\,s^{-2}}\) corresponds to changing speed by  
  \(\Delta u \approx 8.6\times 10^{-3}\,\mathrm{m\,s^{-1}}\) per day.

### 7.2 Stress magnitude relative to wind stress

Take typical values:

- \(\rho_i \approx 900\,\mathrm{kg\,m^{-3}}\)
- \(h \approx 1\,\mathrm{m}\)
- \(a \approx 1\)
- so \(m \approx 900\,\mathrm{kg\,m^{-2}}\)

Then:

\[
m\,C_s \approx 900 \times 10^{-4} = 0.09\,\mathrm{Pa}
\]

For \(F_2\approx 1\), this yields a peak drag scale comparable to typical Antarctic wind stresses (~0.1 Pa). This is consistent with the motivation used by Liu et al. (see Section 11).

### 7.3 A practical “solve for C\_s” tuning rule

If you have a target peak lateral stress magnitude \(\tau_*\) (Pa) for a representative coastal cell:

\[
C_s \approx \frac{\tau_*}{m\,F_2}
\]

Example:

- \(\tau_* = 0.1\,\mathrm{Pa}\)
- \(m = 900\,\mathrm{kg\,m^{-2}}\)
- \(F_2 = 0.5\)

\[
C_s \approx \frac{0.1}{900\times 0.5} \approx 2.2\times 10^{-4}\,\mathrm{m\,s^{-2}}
\]

This is close to the AFIM configuration values currently being tested (order 10\ :sup:`-4`).

---

## 8. Resolution dependence: why **C\_s** can vary by 10–1000×

Users frequently observe that a “reasonable” \(C_s\) at kilometre-scale resolution can be far too weak at coarse climate-model resolutions, and vice versa. There are two main reasons.

### 8.1 Area-averaging a line force

Coastal drag ultimately originates from interactions along (or very near) a boundary—effectively a **line** or narrow **strip** process. Models, however, apply momentum tendencies as **area-averaged** stresses over grid cells.

If the effective coastal boundary layer has physical width \(w\) (km), but the grid spacing is \(\Delta\) (km), the fraction of a coastal cell that “should feel” boundary drag is roughly:

\[
f \sim \frac{w}{\Delta}
\]

If the scheme applies drag over the whole cell, you may need to scale the *area-mean* stress by approximately \(1/f\sim\Delta/w\) to preserve the same **integrated** momentum sink.

Since \(\tau \propto C_s\), a first-order expectation is:

\[
C_s(\Delta) \propto \frac{\Delta}{w}
\]

### 8.2 Sub-grid coastline complexity (partly handled by F2)

The form factor \(F_2\) is designed to represent unresolved coastline roughness and orientation within a cell. When \(F_2\) is computed from high-resolution coastline geometry, it can compensate for some implementation-level resolution effects (coarse cells may contain more sub-grid coastline length, giving larger \(F_2\)).

However, **F2 cannot fully remove resolution dependence**, because:

- the drag still enters as an **area-mean** stress,
- landfast ice physics depends on the ability to form and maintain narrow anchored features,
- and the effective width \(w\) of the coastal shear/anchor zone is not constant across regimes.

### 8.3 Worked examples (order-of-magnitude)

Assume:

- \(w = 5\,\mathrm{km}\) is the effective width of the coastal fast-ice anchor/shear zone.

Then:

| Grid spacing \(\Delta\) | \(\Delta/w\) | Implication for \(C_s\) |
|---:|---:|---|
| 5 km  | 1   | baseline tuning (e.g., 1–3×10\ :sup:`-4` m s\ :sup:`-2`) |
| 50 km | 10  | 10× larger \(C_s\) may be needed if \(F_2\) saturates |
| 500 km| 100 | 100× larger \(C_s\) may be needed (very coarse GCM) |

If instead \(w\sim 0.5\,\mathrm{km}\) (e.g., narrow fjords or grounded-iceberg pinning), then:

- 50 km grid: \(\Delta/w\sim 100\)
- 500 km grid: \(\Delta/w\sim 1000\)

This yields the **10–1000×** scaling range often encountered in practice.

**Important:** in AFIM, because \(F_2\) is computed from high-resolution coastline/GI geometry, some of the above scaling may be moderated. The table should be treated as a **tuning heuristic**, not a strict rule.

---

## 9. Verification checklist

When you enable free-slip + CDP, verify in this order:

1. **Configuration check**
   - `boundary_condition = 'free_slip'` in `&dynamics_nml`
   - `grid_ice = 'C'` (or compatible CD/C-grid path)

2. **Runtime log**
   - confirm `boundary_condition` is read correctly (consider adding a one-line `nu_diag` print in init for traceability)
   - confirm CDP is enabled (`coastal_drag = .true.`)

3. **Form-factor loading**
   - check `nu_diag` output for:
     - F2 file name
     - mapping method
     - min/max of `F2E` and `F2N`
     - nonzero counts (to ensure F2 is not being zeroed everywhere)

4. **Stress diagnostics**
   - confirm `Kux*`/`Kuy*` diagnostics are nonzero in coastal zones
   - confirm `KuE/KuN` are zero over land and nonzero in intended coastal faces

5. **Behavioural sanity checks**
   - compared to no-slip, near-coast velocities should not be artificially pinned *everywhere*
   - with CDP enabled, landfast zones should exhibit a controlled damping consistent with \(C_s\) magnitude

---

## 10. Caveats and known limitations

- The free-slip implementation is applied **within the U-point strain-rate operator** for the C-grid EVP pathway. If other parts of the code still explicitly zero velocities adjacent to land (or deactivate points due to B-grid geometry), those effects can dominate.
- The discrete “reflection/copy” approach is an approximation to free-slip on a **grid-aligned** boundary. For highly oblique coastlines relative to the grid, the effective tangential/normal decomposition is only approximate.
- When combining free-slip with CDP, you should interpret results as a *controlled numerical/physical experiment*: free-slip changes what the solver permits near boundaries; CDP then imposes the intended physical damping.

---

## 11. References

- Liu, Y., Losch, M., Hutter, N., & Ma, L. (2022). *A new parameterization of coastal drag to simulate landfast ice in deep marginal seas in the Arctic.* **JGR: Oceans**, 127. doi: 10.1029/2022JC018413.  
- Losch, M., Menemenlis, D., Campin, J.-M., Heimbach, P., & Hill, C. (2010). *On the formulation of sea-ice models. Part 1: Effects of different solver implementations and grid types.* **Ocean Modelling**, 33, 129–144.  
- Bouillon, S., Fichefet, T., Legat, V., & Madec, G. (2009). *The elastic–viscous–plastic method revisited.* **Ocean Modelling**, 27, 174–186.
