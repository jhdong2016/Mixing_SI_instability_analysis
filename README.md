# Mixing and Symmetric Instability

Companion documentation for the manuscript:

**Effect of Vertical Mixing on the Horizontal Scale of Symmetric Instability in the Ocean Surface Boundary Layer**  
**Author:** Jihai Dong  
**Institution:** Nanjing University of Information Science and Technology, Nanjing, China  
**Contact:** jihai_dong@nuist.edu.cn

## 1. Overview

The study examines how prescribed vertical mixing affects the growth rate and fastest-growing horizontal wavelength of symmetric instability (SI) at a fixed thermal-wind front. The analysis uses nondimensional, nonhydrostatic, linearized Boussinesq equations with constant vertical eddy viscosity and buoyancy diffusivity.

The background shear and stratification remain fixed within each mixing sweep. The baseline contains no horizontal diffusion. The calculations therefore isolate the direct effects of vertical mixing on perturbations, rather than the effects of an evolving background front or a nonlinear feedback between SI and mixing.

## 2. Numerical methods and reference configuration

> **All numerical results in the main text, including Figure 3 and the data used to fit Equation (7), use the 80-cell staggered-grid solver. The 96-interval Chebyshev solver provides an independent robustness check; it is not the source of the main-text Figure 3 results.**

### Staggered-grid solver

Horizontal velocity, buoyancy, and pressure are located at cell centers, while vertical velocity is located at cell interfaces. Pressure is eliminated by projecting the discrete velocity onto the divergence-free subspace. At each horizontal wavenumber pair, the eigenvalue with the largest real part determines the growth rate.

For finite mixing, the **perturbations** satisfy impermeable, stress-free, and insulating boundary conditions:

$$
w=0,\qquad u_z=v_z=b_z=0\qquad\text{at }z=-1,0.
$$

These boundary conditions apply to the perturbations, not to the total prescribed background fields.

### Baseline parameters

| Quantity | Manuscript configuration |
| --- | --- |
| Vertical discretization for main-text numerical results | Staggered grid with 80 cells |
| Background velocity | $U(z)=z+1/2$, $V=W=0$ |
| Background buoyancy | $B(y,z)=z-Ri^{-1}y$ |
| Vertical domain | $-1\le z\le0$ |
| Nonhydrostatic parameter | $\delta=H/L=0.10$ |
| Turbulent Prandtl number | $Pr=1$ |
| Vertical Ekman number | $E_v=\nu_v/(fH^2)$ |
| Buoyancy-diffusion parameter | $E_b=E_v/Pr$ |
| Horizontal diffusion | None |
| Strict-SI branch | $k=0$, $\ell>0$ |
| Richardson numbers in Figure 3 | 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80 |
| Mixing range in Figure 3 | $10^{-5}\le E_v\le10^{-2}$ |

Here, $k$ and $\ell$ are the nondimensional along-front and cross-front wavenumbers, respectively. The velocity and horizontal length scales are $U_0=S_0H$ and $L=U_0/f$, respectively, and the gradient Richardson number is $Ri=N_b^2/S_0^2$.

### Independent numerical checks

Supporting Information Text S1 and Figure S1 compare 80-, 160-, and 320-cell staggered-grid solutions with an independent 96-interval Chebyshev solver. Peak wavenumbers are optimized independently for each resolution and solver.

Across the 40 paired growing cases reported in Text S1, the maximum absolute relative differences between the 80-cell staggered-grid and 96-interval Chebyshev solutions are **0.297% for peak growth rate** and **10.82% for selected wavelength**. These percentages characterize differences between methods for that comparison sample; they are not universal error bounds.

A separate comparison of independently optimized 64- and 96-interval Chebyshev solutions gives a maximum selected-wavelength difference of **0.09% across 78 growing cases**. This is a convergence test of the independent Chebyshev reference solutions, not a change to the main-text solver.

## 3. Growth-rate and wavelength diagnostics

Perturbations have the normal-mode form

$$
\exp\!\left[i(kx+\ell y)+\sigma t\right].
$$

All reported growth rates are nondimensional. Positive $\operatorname{Re}(\sigma)$ denotes exponential growth; the dimensional growth rate is $f\operatorname{Re}(\sigma)$.

For the strict-SI branch, the peak growth rate and selected wavelength are

$$
\sigma_*(E_v)=\max_{\ell>0}\operatorname{Re}\!\left[\sigma(0,\ell;E_v)\right],
\qquad
\frac{\lambda_*}{L}=\frac{2\pi}{\ell_*}.
$$

A selected wavelength is reported only when the maximum has positive growth and lies inside the scanned wavenumber interval. Do not assign a selected wavelength to a non-growing case or interpret a scan-edge maximum as a resolved interior peak. In the inviscid strict-SI limit, the growth rate approaches a finite short-wave limit without selecting a finite fastest-growing wavelength.

## 4. Correspondence with manuscript figures

The figures use the following configurations, as described in the manuscript and Supporting Information.

| Figure | Analysis and configuration |
| --- | --- |
| **Figure 1** | Two-dimensional growth-rate spectra at $Ri=0.25$, 0.70, and 0.95, comparing $E_v=0$ with $E_v=10^{-3}$. Numerical results use the 80-cell staggered grid. |
| **Figure 2** | Strict-SI growth-rate spectra and fastest-growing vertical-velocity structures at $Ri=0.70$. Finite-mixing cases use $E_v=10^{-5}$, $10^{-4}$, $10^{-3}$, and $3\times10^{-3}$ with the 80-cell staggered grid. The dashed inviscid reference curve is analytical. |
| **Figure 3** | Peak strict-SI growth rate and normalized selected wavelength across the eight Richardson numbers and mixing range listed above. The results and the 78 growing cases used for the empirical fit come from the 80-cell staggered-grid solver. |
| **Figure 4** | Energy-budget and velocity-weighted vertical-gradient diagnostics at strict-SI peaks for $Ri=0.70$, using the 80-cell staggered-grid solver. |
| **Figure S1** | Vertical-resolution and solver comparisons: 80, 160, and 320 staggered-grid cells versus a 96-interval Chebyshev reference. The 40-case comparison uses eight Richardson numbers and $E_v=10^{-5}$, $3\times10^{-5}$, $10^{-4}$, $3\times10^{-4}$, and $10^{-3}$. |
| **Figure S2** | Sensitivity to $\delta=0.05$, 0.10, and 0.20 at $Ri=0.70$ and $Pr=1$, using the 80-cell staggered-grid solver. |

Supporting Information Text S1 describes the numerical methods and convergence tests; Text S2 presents the weak-mixing interpretation of wavelength selection; Text S3 defines the energy and vertical-gradient diagnostics.

## 5. Reported wavelength scaling

For the 78 growing cases in the main-text Figure 3 sample, the manuscript reports the empirical approximation

$$
\frac{\lambda_*}{L}\approx3.33\,Ri^{0.83}E_v^{1/4}.
$$

The fit is performed in logarithmic space with the mixing exponent fixed at $1/4$, for $\delta=0.10$ and $Pr=1$. Within the fitted sample, the mean and maximum absolute relative errors are **3.6%** and **24.5%**, respectively, with larger departures under strong mixing near suppression.

This relation summarizes the sampled growing cases. It is not an exact dispersion relation, does not predict the suppression threshold, and does not establish that its coefficients are independent of $\delta$ or $Pr$.

The near-quarter-power dependence is interpreted as a weak-mixing balance between a finite-depth growth deficit proportional to $\ell^{-2}$ and an approximate vertical-mixing attenuation proportional to $E_v\ell^2$. See Text S2 for the assumptions and omitted effects.

## 6. Energy and vertical-gradient diagnostics

The total perturbation-energy budget is

$$
2\operatorname{Re}(\sigma)E=P_S+P_F-\epsilon_v-\epsilon_b,
$$

where $P_S$ is shear production, $P_F$ is potential-energy production associated with the cross-front buoyancy flux and background lateral buoyancy gradient, and $\epsilon_v$ and $\epsilon_b$ are nonnegative viscous and buoyancy-diffusive losses. Vertical buoyancy conversion transfers energy between perturbation kinetic and potential energy and cancels from their sum; it is not the frontal-production term $P_F$.

Each signed budget term in Figure 4a is divided by **twice the total perturbation energy**, $2E$. The velocity-weighted vertical-gradient metric in Figure 4b is defined in main-text Equation (9), with its discrete implementation in Text S3. Its contribution to viscous attenuation satisfies

$$
\frac{\epsilon_v}{2E_K}=E_v m_{\mathrm{eff}}^2,
\qquad
-\frac{\epsilon_v}{2E}=-E_v m_{\mathrm{eff}}^2\frac{E_K}{E}.
$$

The kinetic-energy and total-energy normalizations should not be interchanged.

## 7. Reference values and reproduction notes

The following rounded values are reported in the manuscript for strict SI at $Ri=0.70$, $\delta=0.10$, and $Pr=1$, using the 80-cell staggered-grid solver. They are reference values from the manuscript, not results of a new validation run.

| $E_v$ | Peak growth rate $\sigma_*$ | Selected wavelength $\lambda_*/L$ |
| --- | --- | --- |
| $10^{-4}$ | 0.578 | 0.244 |
| $10^{-3}$ | 0.394 | 0.442 |
| $3\times10^{-3}$ | 0.112 | 0.588 |

Both `main_SI_instability.m` and the implicit method default in `solve_SI_instability.m` select the **80-cell staggered-grid implementation**. The main-text configuration uses `opt.method = 'staggered'`, `opt.N = 80`, and `opt.scan_N = []`, so the coarse wavenumber scan and final peak refinement use the same vertical resolution. Chebyshev calculations remain available as explicitly selected independent checks.

Keep staggered-grid and Chebyshev outputs clearly distinguished, including their solver type, vertical resolution, physical parameters, and peak-search settings. The Figure 3 empirical fit must use the main-text staggered-grid results, not a substituted Chebyshev sample.

The repository includes the seven processed MATLAB data files and six figure-production scripts listed below. Replotting these archived data is distinct from recalculating the eigenvalue problem. Record the Git commit or release used for either task, together with the MATLAB version and run settings.

### Included figure inputs

| Figure script | Required files in `data/` | Output basename in `output/` |
| --- | --- | --- |
| `Figure1.m` | `Figure1_data.mat` | `Figure1_N80` |
| `Figure2.m` | `Figure2_data.mat`, `Figure2_SI_vertical_modes_data.mat` | `Figure2_N80` |
| `Figure3.m` | `Figure3_data.mat` | `Figure3_N80` |
| `Figure4.m` | `Figure4_data.mat` | `Figure4_N80` |
| `FigureS1.m` | `FigureS1_data.mat` | `FigureS1` |
| `FigureS2.m` | `FigureS2_data.mat` | `FigureS2` |

### Replot the archived figures

Set the MATLAB current folder to the repository root and run:

```matlab
Figure1
Figure2
Figure3
Figure4
FigureS1
FigureS2
```

Each script reads the supplied MAT files and saves PDF, TIFF, and MATLAB FIG output in `output/`. These scripts do not solve the eigenvalue problem and do not modify the archived input data. They use MATLAB graphics functions including `exportgraphics` and `tiledlayout`; use a MATLAB release supporting these functions.

### Recalculate the linear stability problem

Run:

```matlab
main_SI_instability
```

The default demonstration computes strict-SI peaks at `Ri = 0.70` and `Ev = [1e-4, 1e-3, 3e-3]`. To recalculate the complete Figure 3 parameter sweep, use the eight-element `Ri_values` and eleven-element `Ev_values` vectors documented near the top of `main_SI_instability.m`, while retaining the staggered-grid settings above.

The optional fixed-wavenumber benchmark block explicitly uses Chebyshev96 for an **independent reference check**. Those benchmark calculations do not change the method used by the main parameter sweep.

Recalculated results are saved in `SI_solver_output/`, including `SI_results.mat`, `SI_peak_summary.csv`, and, for strict-SI runs, `SI_numerical_data.csv`. The solver driver does **not** automatically replace the seven archived `data/*.mat` files or rebuild every figure's data package. The plotting scripts continue to use the archived MAT files until those inputs are explicitly regenerated and replaced.

## 8. Manuscript and archive

The mathematical formulation, interpretation, and figure definitions are documented in the manuscript and its Supporting Information. The Open Research statement specifies that the manuscript-matched repository version will be archived in Zenodo and assigned a DOI before publication. Reference the specific code version used in subsequent analyses.

For questions about the manuscript configuration or reproduction materials, contact Jihai Dong at jihai_dong@nuist.edu.cn.
