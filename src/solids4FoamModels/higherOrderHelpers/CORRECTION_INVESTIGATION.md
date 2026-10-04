# Handoff: quadratic-reconstruction order loss in solids4foam

Date: 2026-09-17. Prepared for independent review, including Claude.
This is an experimental numerical-method investigation, not a validated fix.

## Question to investigate

Why does a curvature-difference correction improve displacement convergence
from approximately second to third order in the standalone 1-D example, but
fail to restore third order reliably in the 3-D solid solver? Please inspect
the implemented operators rather than assuming either the correction formula
or its implementation is correct. In particular, independently verify beta,
boundary treatment, quadrature, derivative recovery and the complete residual.

The user wants the correction separate from alpha stabilisation, added once
per face as extra traction. Do not modify alpha, fGrad or faceCentreValues
silently to improve the observed convergence.

## Repository and reproducibility

- Checkout: `/Volumes/OpenFOAM/work/solids4foam-ib/solids4foam-ib-2412`.
- Branch: `Feature-orderFixThirdOrder`.
- HEAD at handoff: `202bed96d5771a80138c5147e7146853d8007b0f`.
- Important: experimental changes are uncommitted, including the untracked
  `quadraticGradientCorrection` directory. HEAD alone does not reproduce them.
- Runtime tested: OpenFOAM-v2412 on macOS, serial, PETSc SNES, linearElastic.
- Source `/Volumes/OpenFOAM/OpenFOAM-v2412/etc/bashrc`, then explicitly set
  `SOLIDS4FOAM_DIR` to the checkout above; the inherited value can be wrong.
- Build library from `src/solids4FoamModels` using `./Allwmake -j 6`.
- User-owned higherOrderHelpers README and TO_DO_LIST files also exist;
  preserve unrelated edits. No commits were made during these experiments.

## What p means and which unknowns are stored

- p=1: linear reconstruction; nominally second-order displacement accuracy.
- p=2: quadratic reconstruction; desired third-order displacement accuracy.
- p=3: cubic reconstruction; nominally fourth-order displacement accuracy.
- `movingLeastSquares`: unknowns are point values at cell centres.
- `kExactLeastSquares`: unknowns are cell averages. Taylor reconstruction
  subtracts the cell moments so that its average equals the stored unknown.
- Manufactured displacement errors compare point values for MLS and analytical
  cell averages for k-exact. Stress comparisons remain pointwise in these cases.

The reported order problem concerns p=2. Neither the correction nor its
geometry-only cubic probes change the displacement reconstruction to p=3.

## Standalone 1-D reference

Directory: `/Volumes/OpenFOAM/work/run/highOrderOddOrdersAnalysis`.
Read `README.md` and `1DMMS_FV_correction.py` there.
These files were previously under a `papers` subdirectory; older audit scripts
still refer to that location and need their import path adjusted to rerun.

This is a cell-centred finite-volume scalar problem with point-value unknowns,
exact cell source integrals and a shared conservative face flux. The original
face derivative uses a weighted quadratic fit on an eight-cell interior
stencil. Independent seven-cell quadratic fits provide cell second derivatives.

```text
g_original = sum_j c_j * u_j
beta_f = sum_j c_j * (x_j - x_f)^3 / 6
thirdDerivative_f = (secondDerivative_R - secondDerivative_L)/(x_R - x_L)
g_corrected = g_original - beta_f * thirdDerivative_f
```

Boundary faces use the nearest two interior cell second derivatives. The
script also contains a different four-point third-derivative correction;
do not confuse it with the `curvature` variant discussed here.

The documented 1-D curvature-corrected L2 orders are approximately 2.98 on
uniform meshes and 3.10 on one seeded perturbed-mesh family. These are prior
standalone results, not 3-D validation. There is no alpha stabilisation term
in this 1-D script.

## Important 3-D source files

Paths below are relative to `src/solids4FoamModels`:

- `higherOrderHelpers/quadraticGradientCorrection/quadraticGradientCorrection.C`
  and `.H`: experimental geometry coefficients, derivative recovery and traction.
- `higherOrderHelpers/quadraticGradientCorrection/README.md`: detailed equations.
- `higherOrderHelpers/leastSquaresScheme/leastSquaresScheme.C` and `.H`:
  switch, scale and lazy correction-object ownership.
- `higherOrderHelpers/movingLeastSquares` and `kExactLeastSquares`:
  original reconstruction coefficients and quadrature.
- `higherOrderHelpers/fvMeshQuadrature`: face/cell quadrature and moments.
- `solidModels/linGeomTotalDispSolid/linGeomTotalDispSolid.C`:
  `addHighOrderCorrection()` and its call in the PETSc residual.

The focused utility is `applications/test/quadraticGradientCorrection`.
The separate audit executable described below lives outside the repository.

## Current correction configuration and scope

Inside `linearGeometryTotalDisplacementCoeffs.highOrderCoeffs.displacement`:

```text
type                     kExactLeastSquares; // or movingLeastSquares
polynomialOrder          2;
curvatureCorrection      true;
curvatureCorrectionScale 1;                  // theta, range [0,1]
```

The correction switch defaults to false. Zero theta bypasses the correction
object entirely and reproduces the original residual. Current active scope:
serial 3-D, one constant-property linearElastic law, PETScSNES,
highOrderJacobian=false, solvePressure=false, no coupled or symmetry patches.
Parallel, nonlinear material and analytical-Jacobian correction support have
not been implemented.

## Current beta construction

For each original face quadrature point q, apply the original gradient
coefficients to the ten factorial-scaled cubic coordinate terms centred on q:

```text
rx^3/6, rx^2*ry/2, rx^2*rz/2, rx*ry^2/2, rx*ry*rz,
rx*rz^2/2, ry^3/6, ry^2*rz/2, ry*rz^2/2, rz^3/6
```

Their exact gradients at q are zero, so the reconstructed gradients give
the corresponding cubic-error responses. Each beta is a vector. MLS samples
point values; k-exact samples cell averages using geometric moments.
K-exact boundary-data coefficient contributions are included in beta.

```text
faceBeta[m] = sum_q faceWeight[q] * beta[q][m] / faceArea
deltaGradient(u)_f = -sum_m faceBeta[m] * estimatedThirdDerivative_f[m]
```

The helper owns auxiliary cell quadrature of degree p+1 for geometric moments.
It does not raise the original solver quadrature. The original reconstruction
requests face quadrature order p-1, hence order 1 for p=2. The existing
polynomial traction test compares with the analytical traction evaluated using
that same rule, not independently exact higher-order face integration.

## Current p=2 derivative recovery: direct two-cell variant

```text
d = |x_N - x_P|
e = (x_N - x_P)/d
J = (secondDerivative_N - secondDerivative_P)/d
```

The second derivatives come from the existing cell reconstruction coefficients.
Their calculation still uses each cell's stencil, but J uses only the two cells.
In a local frame with x along e, J determines xxx, xxy, xxz, xyy, xyz and xzz.
It cannot determine yyy, yyz, yzz or zzz. Those contributions are not corrected.

The code constructs the directional tensor in global coordinates without a
local coordinate transformation. With v=J*e and s=e*v:

```text
thirdDerivative_ijk = e_i*J_jk + e_j*J_ik + e_k*J_ij
                    - e_i*e_j*v_k - e_i*e_k*v_j - e_j*e_k*v_i
                    + e_i*e_j*e_k*s
```

Contracting the last index with e returns J. Purely transverse components are
zero; reversing the cell pair leaves this completed tensor unchanged.
There is no p=2 stencil-wide fit or response calibration in the current code.

At a physical boundary, choose the owner and its face-connected neighbour
most aligned with the inward face normal. If none lies inward, choose the
nearest inward centre from the owner's existing reconstruction stencil.
This fallback still uses only two cell derivatives. This one-sided treatment
is experimental and has no established boundary-order guarantee.

The p=1 correction was left unchanged: a stencil-wide recovery of second
derivatives from first derivatives, with a six-by-six response calibration.

## Traction and solver integration

Assemble the three displacement-component gradient increments into deltaGradD:

```text
deltaSigma = 2*mu*symm(deltaGradD) + lambda*tr(deltaGradD)*I
deltaTraction = theta * (normal & deltaSigma)
```

Mu and lambda come from the actual linearElastic law. Stiffness is already
included; another impK multiplier would count it twice. The correction returns
traction, not force. The existing divergence applies face area once.

Residual sequence:

1. Original gradients, constitutive law and quadrature-integrated traction.
2. Original alpha stabilisation traction.
3. Separate high-order correction traction.
4. Enforce prescribed traction boundary conditions.
5. Assemble face-force divergence and source/time terms.

For constant material coefficients and one normal per face, applying a common
derivative estimate with point-dependent beta at every original quadrature point
is equivalent to using the weighted face-average beta once per face. Moving the
same addition into the quadrature loop alone would not change the residual.

## Previous variants: do not mix their results

Earlier experimental versions also modified fGrad, faceCentreValues and/or
quadrature. Their improved orders cannot be attributed to the current separate
traction correction. Those hidden corrections have been removed.

A later face-only variant recovered all ten third derivatives from stencil-wide
changes in second derivatives, then applied an inverse ten-by-ten geometry
response matrix. This reduced errors but did not uniformly restore order.
Its source snapshot is in the study directory below. It is not the current p=2
implementation. The current implementation replaced only this recovery step
with the direct two-cell method described above.

## MMS studies before the alpha-zero experiment

Parent directory: `/Volumes/OpenFOAM/work/run/MMS_orderTesting`.

- `curvature-convergence.VIzn1v`: original references; p=2 uncorrected cases
  have `false` in their names. Avoid the older `true` cases when comparing.
- `face-traction-mms.khCc1g`: face-only correction with stencil-wide recovery.
- `directional-mms.2Nwtzs`: current two-cell correction, alpha=0.1, theta=1.

Each includes logs and/or analysis files; the latter two have README reports,
plots, CSVs and source snapshots. Six meshes per family:

```text
hex: 125, 216, 343, 512, 1000, 1728 cells
tet: 703, 1561, 2577, 3510, 6290, 9963 cells
```

Domain dimensions are 0.2 in each direction. Slopes use h=Ncells^(-1/3),
equivalent up to a constant to a domain-scaled characteristic length.
The table gives displacement L2 slopes fitted on the finest three meshes:

| Method | Mesh | Original | Full stencil recovery | Two-cell recovery |
| --- | --- | ---: | ---: | ---: |
| MLS | Hex | 2.836 | 3.225 | 1.966 |
| k-exact | Hex | 2.813 | 2.684 | 1.983 |
| MLS | Tet | 2.621 | 2.667 | 2.584 |
| k-exact | Tet | 2.606 | 2.432 | 2.502 |

All 24 two-cell runs converged at alpha=0.1. Finest hex errors increased relative
to the original solver. Zero theta controls reproduced saved displacement
values and printed displacement/stress errors exactly. A successful linear or
nonlinear solve is not evidence of the desired discretisation order.

## Cubic diagnostic findings

Directory: `/Volumes/OpenFOAM/work/run/MMS_orderTesting/audit-directional.5xpxU0`.
Read its README, `analysis.txt`, `log.1d`, `log.final.hex6`, `log.final.hex3`.
Source and scripts are included. No production code was changed by this audit.

For u=x^3, exact u_xx=6*x and u_xxx=6. On the central line, first pair estimates
are approximately 1.590 in the actual 1-D script, 0.270 in 3-D MLS and 0.224 in
3-D k-exact. Interior pairs recover 6. Boundary-layer pair errors persist with
refinement. The measured weights explain the cell derivative errors:

```text
E_P = reconstructed_u_xx_P - 6*x_P
J_xx = 6 + (E_N - E_P)/(x_N - x_P)
```

The 1-D code also has boundary-biased derivatives despite its improved global
convergence. Thus this weakness alone does not explain the 1-D/3-D difference.

Independent beta evaluation with analytical cell averages gave these finest-hex
traction errors for the synthetic vector field (1,-0.7,1.3)*x^3:

| Exact-input diagnostic | MLS | k-exact |
| --- | ---: | ---: |
| Exact second derivatives, x-aligned faces | 6.7e-16 | 3.9e-15 |
| Exact second derivatives, all internal faces | 1.82e-3 | 1.67e-3 |
| Full exact third derivative, all internal faces | 1.0e-15 | 5.6e-15 |

For a y-directed pair, exact second derivatives of x^3 are identical, so its
directional correction is zero. Its face-gradient stencil can still have an
x^3 error, and elastic traction depends on the full gradient. This establishes
a limitation of the two-cell completion. It does not establish that those
individual face errors survive cancellation in the complete cell residual.

The exact-derivative substitution covered internal faces for x^3, not all
cubic modes or physical boundary faces. No simple internal sign, stiffness or
factorial error was found for that mode; general correctness remains unproven.

## Validation limits

- v2412 build succeeds. Other OpenFOAM versions were not tested in this study.
- Forty-eight polynomial method/order/mesh checks pass, but for current p=2
  these assert degree-2 exactness, scale linearity and unchanged original
  operators. Cubic errors are diagnostics, not passing cubic-exactness checks.
- A separate rank-one cubic identity validates directional tensor projection
  and pair reversal, including rotated directions, to about 9e-16.
- The test uses the original face quadrature for analytical reference traction.
- No full-residual polynomial decomposition or complete beta mapping audit has
  yet been implemented. Do not describe these proposed tests as completed.

## Alpha-zero experiment

Directory: `/Volumes/OpenFOAM/work/run/MMS_orderTesting/alpha-zero-mms.6vu6z5`.
The new cases copy `directional-mms.2Nwtzs` and change only
`linearGeometryTotalDisplacementCoeffs.stabilisation.momentum.scaleFactor`
from 0.1 to 0. The correction remains enabled, theta=1, p=2. Original meshes,
case sources, Allrun scripts and solver implementation were not modified.

The same PETSc settings were retained: matrix-free SNES, LGMRES/Hypre,
KSP limit 1000, SNES limit 30, snes_rtol=1e-10 and snes_stol=1e-12.
Removing alpha means its computed traction is multiplied by zero; the alpha
class and original reconstruction path remain selected.

| Method | Mesh | Order, alpha=0.1 | Order, alpha=0 |
| --- | --- | ---: | ---: |
| MLS | Hex | 1.966 | 2.427 |
| k-exact | Hex | 1.983 | 2.430 |
| MLS | Tet | 2.584 | 3.161 |
| k-exact | Tet | 2.502 | No converged solves |

All orders use the finest three meshes, with the two-cell correction enabled
on both sides of the comparison. Eighteen cases converged. All six k-exact
tet cases hit DIVERGED_ITS at 1000 KSP iterations, followed by
DIVERGED_LINEAR_SOLVE at zero nonlinear iterations. No errors or slopes are
assigned to those failed solves. This is failure under the tested solver and
iteration limit, not proof of singularity or impossibility with another solver.

For alpha=0, the finest-pair slopes were 2.313 (MLS hex), 2.349 (k-exact hex)
and 2.892 (MLS tet). Thus the fitted 3.161 on tet is not a uniform pairwise
rate or a demonstrated asymptotic order. Finest displacement L2 errors:

| Method | Mesh | Alpha=0.1 | Alpha=0 |
| --- | --- | ---: | ---: |
| MLS | Hex | 1.54421e-9 | 1.48171e-9 |
| k-exact | Hex | 1.72761e-9 | 1.61848e-9 |
| MLS | Tet | 4.38231e-10 | 4.51930e-10 |

Alpha influences both apparent order and solve robustness, but removing it
does not restore third order consistently. A larger fitted slope can coexist
with a larger error on the finest mesh (MLS tet). An alpha=0, correction-off
control was not run in this experiment, so these data alone do not isolate
the correction's benefit in the absence of alpha.

The study includes `convergence.png`, `convergence.csv`, `slopes.txt`, all
case logs, a source snapshot and scripts to reproduce the run and analysis.

## Interior-vs-boundary decomposition experiment

Directory: `/Volumes/OpenFOAM/work/run/MMS_orderTesting/interior-boundary-split.Ne38yG`.
Tested whether the order loss is a coarse-mesh boundary-shell artifact: these
3D meshes have 42-78% boundary-adjacent cells at the sizes tested (5-12 cells
per side), versus 0.1-6% boundary cells at the 1D script's tested sizes
(32-2048 cells), so a boundary-localised deficiency could plausibly dominate
the 3D fitted order without affecting the 1D result. No new solves were run;
this reprocesses already-saved `1/DDifference` fields from
`curvature-convergence.VIzn1v`, `directional-mms.2Nwtzs` and
`alpha-zero-mms.6vu6z5`, split into boundary-adjacent cells (any cell owning
a boundary face) and interior-only cells.

The metric is unweighted per-cell RMS, not the solver's volume-weighted L2
norm, so this is a diagnostic approximation, not a recomputation of the
logged order; a sanity check against the logged k-exact hex correction order
(1.983) reproduced 1.9835 by this method, so the pipeline is trusted for
comparing all-cells versus interior-only on the same footing.

| Study | Method | Mesh | Order, all cells | Order, interior-only |
| --- | --- | --- | ---: | ---: |
| Baseline (uncorrected) | MLS | Hex | 2.836 | 2.492 |
| Baseline (uncorrected) | k-exact | Hex | 2.813 | 2.262 |
| Baseline (uncorrected) | MLS | Tet | 2.621 | 2.690 |
| Baseline (uncorrected) | k-exact | Tet | 2.606 | 2.657 |
| Correction, alpha=0.1 | MLS | Hex | 1.966 | 2.057 |
| Correction, alpha=0.1 | k-exact | Hex | 1.984 | 2.089 |
| Correction, alpha=0.1 | MLS | Tet | 2.584 | 2.672 |
| Correction, alpha=0.1 | k-exact | Tet | 2.502 | 2.531 |
| Correction, alpha=0 | MLS | Hex | 2.427 | 2.092 |
| Correction, alpha=0 | k-exact | Hex | 2.430 | 2.169 |
| Correction, alpha=0 | MLS | Tet | 3.161 | 3.205 |

Excluding the boundary shell did not move the correction cases meaningfully
closer to third order: hex changed by +0.09 to +0.18, tet by +0.03 to +0.09,
still 2.0-2.7 in every case. For the uncorrected baseline hex cases,
interior-only order was lower than all-cells order (2.836 to 2.492, 2.813 to
2.262), the opposite of what the hypothesis predicts. This is evidence
against the coarse-mesh boundary-shell-fraction explanation as tested here;
the order deficit is not concentrated in the boundary cells, and removing
them does not reveal a hidden third-order interior rate. It does not rule out
boundary treatment as a problem by some other mechanism (for instance,
boundary cells anchored to imposed BC data reading artificially low error
and masking a lower interior order in the all-cells norm), which was not
checked here.

## Quadrature-order experiment

Tests whether under-resolved face quadrature explains the p=2 order loss.
The base reconstruction and the correction's `faceBeta_` construction share
one `fvMeshQuadrature` object (`movingLeastSquares.C`/`kExactLeastSquares.C`,
`makeQuadrature()`); its face order was temporarily raised from
`polynomialOrder_ - 1` (order 1) to `polynomialOrder_ + 3` (order 5, well
past cubic-exact), rebuilt, rerun on the same case matrix as
`directional-mms.2Nwtzs`, then reverted and rebuilt back to the original
state. This is a global change affecting both `fGrad`/`faceCentreValues` and
the correction, not an isolated change to the correction alone. All 24
corrected runs and 4 scale=0 controls converged.

Finest-three fitted order (h = Ncells^-1/3). The baseline uses order-1
quadrature; the other columns give the correction's face quadrature order:

| Method | Mesh | Baseline | Correction, order 1 | Correction, order 5 |
| --- | --- | ---: | ---: | ---: |
| MLS | Hex | 2.836 | 1.966 | 1.952 |
| k-exact | Hex | 2.813 | 1.983 | 1.844 |
| MLS | Tet | 2.621 | 2.584 | 2.521 |
| k-exact | Tet | 2.606 | 2.502 | 2.371 |

Raising quadrature order did not improve any of the four combinations: MLS
hex was essentially flat (-0.014), and k-exact hex (-0.139), MLS tet (-0.063)
and k-exact tet (-0.131) got measurably worse. The scale=0 controls also
changed slightly relative to the original order-1 baseline (e.g. level-1 MLS
hex displacement L2: 2.05911e-08 originally versus 2.08204e-08 at order 5,
about 1.1% different), confirming the base scheme's own face integral was not
fully converged at order 1 either, though the effect is small at this mesh
size. This is evidence against face-quadrature under-resolution as the cause
of the p=2 order loss; raising it made results the same or worse in every
tested case, not better. Full data, the run/analysis scripts and the exact
reverted diff are in
`/Volumes/OpenFOAM/work/run/MMS_orderTesting/quadrature-order-mms.q5rT2k/`.

## Cubic-consistency audit

Implements checks 1 and 2 without solving, in
`/Volumes/OpenFOAM/work/run/MMS_orderTesting/cubic-consistency-audit.I1sEzZ/`.
The actual schemes are applied to 30 cubic displacement modes on all 12 MMS
meshes. Each face traction term is compared with the exact face-integrated
traction: reconstruction `R`, order-1 face quadrature `Q`, alpha term `A` and
the current correction `C`. At p=2, a term that is not exact for cubics is an
O(h^2) face-flux error.

- The correction recomputed from independently measured beta matches
  `addTraction` to at most 2.2e-10 relative. No beta or traction bug was found.
- `R`, `Q` and `A` all scale with slope 2.0-2.3. On the finest meshes, for
  internal faces, `Q/R` is 0.15 (hex) and 0.47-0.51 (tet), and `A/R` is
  0.12-0.16 at alpha=0.1.
- The corrected residual `R+C` is 0.53-0.77 of `R`. With exact cell second
  derivatives in the same pair, 0.51-0.71 of `R` remains.
- Check 2: the owner/other data have rank 6 on every face. Median
  irreducible fractions are 0.68-0.79 (gradient) and 0.48-0.83 (traction).
  Adding the face neighbours of both cells gives rank 10 on every face,
  including boundary faces, with residual below 3.3e-15. The median condition
  number is 3-8 and the maximum is 79.
- 1-D check: adding the solids4foam alpha jump to the reference script gives
  uniform-mesh orders of 1.99 at alpha=0.1 and 2.08 at alpha=0.01, against
  2.98 at alpha=0. Subtracting the cubic response of the one-sided face
  values restores 2.98 at alpha=0.1.

Face-level O(h^2) terms are candidates, not proof, of solution-level order
loss. The face-patch study below found no effect of `Q` on six coarse hex
meshes; the eleven-mesh study further below shows that `Q` does limit the
corrected scheme once the meshes are fine enough. `R+C` (pair limit), `A`
(alpha, shown in 1-D) and `Q` are all real limits.

## Face-patch correction on structured hex, alpha=0

Study: `/Volumes/OpenFOAM/work/run/MMS_orderTesting/facePatch-hex-alpha0.6ebpl4/`.
`curvatureCorrectionRecovery facePatch` fits all ten third derivatives from
the second derivatives of the pair and the face neighbours of both cells,
using their geometry-only cubic response (see the correction README). The
audit and unit test give a corrected cubic traction error at rounding level
for MLS and k-exact, on hex and tet, internal and boundary faces. The
original `twoCell` recovery remains the default and is unchanged.

The hex cases of `directional-mms.2Nwtzs` were rerun with alpha=0 and
`facePatch`, at the original face order 1 and a temporary face order 2
(reverted, library rebuilt):

| k-exact, hex | Finest L2 | Pairwise displacement orders |
| --- | ---: | --- |
| Original, alpha 0.1 | 1.118e-09 | 3.41 3.27 3.12 2.91 2.68 |
| twoCell, alpha 0 | 1.618e-09 | 2.99 2.77 2.64 2.49 2.35 |
| facePatch, alpha 0, order 1 | 1.804e-10 | 4.44 4.69 4.56 4.56 4.40 |
| facePatch, alpha 0, order 2 | 1.260e-10 | 4.17 4.51 4.35 4.41 4.43 |

k-exact reaches at least third order over these meshes at either face order.
Stress converges at second order. MLS with `facePatch` and alpha=0 needed
thousands of linear iterations per Newton step from mesh 2 onwards and did
not converge on meshes 3-6, with a 20000-iteration limit. Its level-2 error
exceeded level 1. The cause of this poor conditioning is not yet known. Six
coarse meshes do not establish the asymptotic order.

With alpha 0.1, `facePatch` also removes the cubic part of the alpha jump
inside its separate traction; `alphaStab` is unchanged (correction README).
The audit confirms that R+C+A is at rounding level on hex and tet faces.
Results use a pinned MMS library; see the study README. Another checkout's
tests overwrote the shared `libmanufacturedSolution` during this study, and
the affected runs were discarded.

| Hex, alpha 0.1 | Finest L2 | Pairwise displacement orders |
| --- | ---: | --- |
| k-exact, uncompensated | 4.317e-10 | 4.43 3.87 3.06 2.32 2.00 |
| k-exact, compensated | 1.311e-10 | 4.46 4.59 4.50 4.46 4.14 |
| MLS, uncompensated | 3.832e-10 | 4.62 3.96 3.47 2.85 2.40 |
| MLS, compensated | 1.544e-10 | 4.69 4.40 4.41 4.38 4.06 |

Uncompensated alpha returns both methods to second order, as in 1-D. With
compensation, MLS also converges robustly (at most 28 linear iterations per
Newton step). Tet meshes have not been run with `facePatch`.

The observed order is close to four, not three, and the displacement is
close to p=3 results. This follows from the design: `facePatch` makes the
face traction exact for cubics, which is p=3's flux consistency. The 1-D
reference reaches only three because of its boundary closure. Using the
exact u''' on the four faces next to each boundary, where its stencils are
shifted, gives uniform-mesh pairwise orders of 4.44, 4.08, 3.93, 3.90, 3.93
and 3.96. The finest error falls from 5.6e-9 to 3.6e-11. One face per
boundary is not enough; that still gives 3.00. The script is
`cubic-consistency-audit.I1sEzZ/oneD/boundary_closure_1d.py`. Removing the
O(h^2) flux error of p=2 requires consistent estimates of all ten third
derivatives. Whether corrected p=2 is worthwhile is therefore a question of
cost against p=3 at equal error.

## Cell-fit correction and tetrahedral meshes

Study: `/Volumes/OpenFOAM/work/run/MMS_orderTesting/p2recovery-cellFit-tet.zzCksW/`.
The default p=2 stencils (64-80 cells, chosen for conditioning) are only
about 12% smaller than p=3's (74-83), and a cubic fit on them is full rank
on every cell of every mesh (`cubic-consistency-audit.I1sEzZ/stencilCheck`).
`curvatureCorrectionRecovery cellFit` therefore fits a cubic to each cell's
own p=2 stencil and averages owner and neighbour at the face, keeping the
residual on the p=2 footprint. It is cubic-exact to rounding with the alpha
compensation, on hex and tet, MLS and k-exact. Alpha 0.1 throughout.

Fitted order / finest-mesh displacement L2:

| Mesh, method | Uncorrected p=2 | cellFit p=2 | p=3 |
| --- | ---: | ---: | ---: |
| Hex, k-exact | 2.81 / 1.12e-9 | 4.53 / 1.34e-10 | 3.48 / 8.4e-11 |
| Hex, MLS | 2.84 / 1.35e-9 | 4.49 / 1.63e-10 | 3.97 / 1.2e-10 |
| Tet, k-exact | 2.61 / 4.47e-10 | 2.84 / 1.55e-10 | 4.10 / 1.5e-11 |
| Tet, MLS | 2.62 / 4.49e-10 | 2.92 / 1.86e-10 | 4.12 / 1.5e-11 |

Fitted slopes use the finest three meshes; tet pairwise orders scatter from
2.3 to 3.6. Stress stays at 1.9-2.1 for every p=2 variant. On tets the
corrected scheme behaves as a third-order scheme, with errors 2.4-2.9 times
below uncorrected p=2, while p=3 on the same meshes shows four. The order
above four is a structured-hex effect. On tets, p=3 is 10-12 times more
accurate than corrected p=2 at 1.2-1.6 times the run time on these small
meshes. `facePatch` on tets is poorly conditioned (about 900 linear
iterations per Newton step for k-exact, failing at the 1000 limit on five
of six meshes); `cellFit` needs at most 31.

## Eleven-mesh study: the face quadrature limits the corrected scheme

Study: `/Volumes/OpenFOAM/work/run/MMS_orderTesting/extended-11levels/`
(README there). Meshes 1 to 11 of the study `Allrun` (hex up to 29,791
cells, tet up to 148,567), alpha 0.1, MLS p=1/2/3 and p=2 with `cellFit`,
k-exact p=1/2/3; every case runs against a private copy of this checkout's
library. Six coarse meshes were not enough to see the asymptotic order.

| MLS, finest three of eleven | Hex | Tet |
| --- | ---: | ---: |
| p=2 | 2.26 | 2.15 |
| p=2 + cellFit, face quadrature order 1 (`p-1`) | 1.90 | 2.51 |
| p=2 + cellFit, face quadrature order 2 | 4.04 | 4.02 |
| p=3 | 4.03 | 3.94 |

With the default face rule the corrected hex order falls from 4.6 on the
coarse meshes to 1.85 on the finest pair, and the corrected error settles
at 9% of the uncorrected p=2 error from mesh 8 on. With a quadratic-exact
face rule the order stays at four on both mesh families; the finest errors
are 1.24 (hex) and 1.68 (tet) times those of p=3 at 75% and 96% of the p=3
run time (cases timed one at a time; the earlier batch timings were
inflated by concurrency). Stress stays second order for every p=2 variant.

This settles the quadrature question: the corrected flux is exact for
cubics at the quadrature points, but the one-point-per-triangle rule
integrates the quadratic traction of a cubic field with an O(h^2) error
(`Q` in the audit). On coarse meshes the correction's own higher-order
remainder hid it, which is why raising the rule earlier changed nothing.
Corrected p=2 therefore needs `faceQuadratureOrder 2`, now a dictionary
entry of the reconstruction. Order 3 gives the same orders and errors
within 10% at higher cost, so order 2 is sufficient. On tets the order-1
rule samples the face centre, so evaluating the correction at the face
centre with one point per triangle is already covered by those runs: it
stalls at 2.5. Evaluating beta once at the face centre instead
(`curvatureCorrectionBeta faceCentre`, MLS only, via the new
`movingLeastSquares::faceGradCoeffsAtPoint()`) gives 1.94 at face order 1
and 2.62 at face order 2 on hex, against 1.90 and 4.04 with the default
quadrature average: the cubic cancellation is exact only when beta is
sampled with the rule that integrates the traction, because the
least-squares weights depend on the evaluation point. k-exact with the
correction has not been run on eleven meshes.

## Q term: closed-form face-integration correction

`curvatureCorrectionFaceIntegration true` keeps one point per triangle and
adds the rule's integration error for a cubic field to the correction:
`0.5*(J_exact - J_rule)/A` contracted with the ten third derivatives from
`cellFit`, where `J` are the second area moments of the face about its
centre, exact and as sampled by the actual quadrature points. The audit
shows `R+C+Q` and `R+Ca+Q+A` at about 1e-17 with the term on. Eleven
meshes, MLS, alpha 0.1, finest-three orders: 4.09 (hex) and 4.05 (tet),
against 1.90 and 2.51 without it and 4.04 and 4.02 with face order 2;
stress 2.0. Finest errors are 1.16 (hex) and 1.21 (tet) times the MLS p=3
errors, at 42% and 77% of the MLS p=3 time (23 s against 55 s on 29,791
hex cells; 97 s against 125 s on 148,567 tets; all timed one at a time).
The Q term costs 5% and 14% over the uncorrected-quadrature run and
replaces the quadratic face rule. k-exact p=3 runs in 23 s and 109 s on the
same meshes with lower errors and third-order stress, so it remains the
competitor to beat. This is the finite-volume form of the flux-quadrature
corrections in Nishikawa, AIAA 2025-3674, eqs. 4.7 and 4.13, with the
third derivatives from the cell fit in place of nodal flux gradients.

## Proposed independent checks

1. Apply all ten cubic monomials to the actual reconstruction. Independently
   compare stored beta with measured fGrad error, including physical boundaries.
2. At a representative face, collect a 6-by-10 matrix J of the actual two-cell
   derivative differences and a 3-by-10 matrix E of gradient errors. Determine
   whether any 3-by-6 mapping K can satisfy K*J=E. Nonzero irreducible residual
   means no coefficients based only on that pair can cancel all cubic modes.
3. Compare original traction, corrected traction, alpha and integrated source
   in the complete residual on exact polynomial displacement, without solving.
   Separate interior and boundary-adjacent cells. Check cancellation across faces.
4. Compare analytical face integration with independently higher-order
   quadrature. A corrected gradient may contain quadratic terms even though
   the original p=2 face rule is requested with order 1.
5. Only after isolating the leading residual error, consider modifying beta,
   recovery, boundaries, quadrature or alpha. Keep changes independently testable.

Open question: a cell-pair correction may be insufficient, but the earlier
full-derivative correction also did not restore order uniformly. The missing
piece may be elsewhere in the force balance, or there may still be an
implementation issue not covered by the single-mode audit.
