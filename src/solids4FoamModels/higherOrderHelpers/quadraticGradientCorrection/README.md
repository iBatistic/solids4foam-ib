# Experimental face-traction correction

This helper supplies a **separate traction contribution once per face**.
It does not modify `fGrad()`, `faceCentreValues()`, alpha stabilisation,
or the quadrature settings used by the original solver. With `facePatch` the
separate contribution also removes the cubic part of the alpha term.

The implementation supports p=1 and p=2. The class retains its original
`quadraticGradientCorrection` name, but its earlier gradient/value-correction
interface has been removed.

## Configuration

In the existing high-order displacement reconstruction dictionary:

```text
type                     kExactLeastSquares;  // or movingLeastSquares
polynomialOrder          2;                   // 1 or 2
curvatureCorrection      true;
curvatureCorrectionScale 1.0;                 // 0 to 1
curvatureCorrectionRecovery twoCell;          // facePatch or cellFit (p=2)
curvatureCorrectionFaceIntegration true;      // Q term, cellFit only
faceQuadratureOrder      1;                   // default p-1; 2 if no Q term
curvatureCorrectionBeta  faceAverage;         // or faceCentre (MLS, diagnostic)
```

`curvatureCorrectionFaceIntegration` adds the face rule's integration error
for cubic fields to the correction (the Q term below), so the default
one-point-per-triangle rule is sufficient; without it, corrected p=2 needs
`faceQuadratureOrder 2`.

`curvatureCorrectionBeta faceCentre` evaluates beta with one reconstruction
of the face stencil at the face centre instead of averaging the quadrature
points. It is a diagnostic: the cubic cancellation is exact only with the
default average, and on eleven hex meshes it gave 1.94 (face order 1) and
2.62 (face order 2) against 1.90 and 4.04 for the default.

The switch defaults to false; the scale defaults to 1 and must be finite
and within [0, 1]. Scaling does not change the geometric beta coefficients.
The recovery defaults to `twoCell`; `facePatch` and `cellFit` are described
below. `faceQuadratureOrder` sets the face rule of the reconstruction
itself (see the results below); it is not changed by the correction.

Current solver integration is in the PETSc residual of
`linearGeometryTotalDisplacement`. A nonzero correction requires one
constant-property `linearElastic` material, `solvePressure false`, a serial
3-D mesh without coupled or symmetry patches, and
`highOrderJacobian false` with a matrix-free Jacobian.

A false switch or zero scale bypasses correction construction and addition.
Original quadrature and alpha remain unchanged, so zero scale recovers the
original residual. The usual restrictions of that original solver still
apply. A nonzero correction has not been implemented for the analytical
Jacobian, parallel runs, symmetry, nonlinear materials or mixed D-p solves.

## Solver sequence

The residual follows this sequence:

1. Evaluate original cell and face-quadrature gradients.
2. Evaluate the original constitutive law and integrate face traction.
3. Add the original alpha stabilisation contribution.
4. Call `addHighOrderCorrection(traction, D)` once to add the face correction.
5. Enforce prescribed traction boundary conditions.
6. Assemble the divergence of face traction multiplied by face area.

```text
total traction = original traction + original alpha traction
               + scale * correction traction
```

The correction uses the displacement passed into the residual, not an
independently cached solution. Internal faces have one shared traction,
giving opposite force contributions to owner and neighbour. Prescribed
traction conditions are applied afterwards and therefore take precedence.
Fixed-displacement faces retain the calculated additional traction.

This design does not add the correction to stored quadrature gradients or
stresses. Stress post-processing still uses the original reconstruction.
No face area is applied twice: the helper returns traction, while the
solver divergence applies area.

## Recovering the missing derivatives

For p=1, recover the six independent second derivatives from spatial changes
in the original reconstructed first derivatives:

```text
raw u_xx = d_x(reconstructed u_x)
raw u_xy = [d_x(reconstructed u_y) + d_y(reconstructed u_x)] / 2
```

For p=1, the spatial changes use a linear, distance-weighted least-squares
fit on the existing cell stencil. The geometry-only response matrix described
below is used only for p=1.

### Direct two-cell correction for p=2

For p=2, calculate the original second derivatives in every cell using
`cellSecondGradCoeffs()`. Then use only the owner P and neighbour N at each
internal face. No additional stencil-wide fit or response calibration is used:

```text
d = |x_N - x_P|
e = (x_N - x_P)/d
J = [secondDerivative_N - secondDerivative_P]/d
```

The second derivatives are constant in each quadratic reconstruction, so J
is shared by all quadrature points on the face. Each cell's second derivative
still uses its original reconstruction stencil, just as in the 1-D example.

With local x along e, the six available third derivatives are:

```text
u_xxx = J_xx;  u_xxy = J_xy;  u_xxz = J_xz
u_xyy = J_yy;  u_xyz = J_yz;  u_xzz = J_zz
```

The pair cannot determine u_yyy, u_yyz, u_yzz or u_zzz. Their contributions
are left uncorrected. This is not complete cubic reconstruction.

The code avoids constructing a local coordinate frame. It forms a symmetric
directional third-derivative tensor in global coordinates. Writing its entries
as `thirdDerivative_ijk`, with i, j and k each denoting x, y or z:

```text
v = J & e
s = e & v
thirdDerivative_ijk = e_i*J_jk + e_j*J_ik + e_k*J_ij
                    - e_i*e_j*v_k - e_i*e_k*v_j - e_j*e_k*v_i
                    + e_i*e_j*e_k*s
```

Contracting its last index with e returns J. Its components with all three
directions perpendicular to e vanish. Thus it implements exactly the six
local-coordinate entries above, not ten independent derivative estimates.
Reversing the pair changes the signs of both e and J, leaving the resulting
third derivative and correction unchanged.

For a physical boundary face, P is the owner. N is its face-connected neighbour
most aligned with the inward face normal.
Selection maximises the positive dot product of the inward area vector with the
unit vector from P to N; exact ties retain the first cell in the adjacency list.
If no immediate neighbour lies inward, use the nearest inward cell centre
from the owner's existing reconstruction stencil. This fallback still selects
only one other cell, without fitting derivative differences over the stencil.
The same two-cell formula is used without extrapolation to the boundary face.
Construction fails if there is no valid inward pair. This is an experimental
one-sided counterpart of using the nearest two cells at a 1-D boundary, not
an assertion of boundary-order accuracy on skewed meshes. Prescribed traction
still overrides the added correction afterwards.

For exact cubic derivatives, J supplies the six directional combinations
exactly. With reconstructed second derivatives there is an additional error:

```text
J_reconstructed = J_exact + (error_N - error_P)/d
```

This error need not cancel on irregular stencils. Also, for a general smooth
field the difference is associated with the cell-pair segment, not necessarily
the face centre. Neither full cubic exactness nor increased solution order is
guaranteed. In 1-D the expression reduces to the original correction:

```text
deltaGradient = -scale * beta * (u_xx_N - u_xx_P)/(x_N - x_P)
```

### Face-patch correction for p=2

`curvatureCorrectionRecovery facePatch` estimates all ten third derivatives.
The patch S of a face is the two-cell pair above plus the face neighbours of
both cells. For a cubic field with third derivatives T, the reconstructed
second derivative of every cell is linear in T. Their deviations from the
patch mean therefore satisfy:

```text
secondDerivative_k - mean_S(secondDerivative) = M_k * T,   k in S
```

M is the 6|S|-by-10 geometry-only response. It is built by applying the
existing second-derivative coefficients to the ten cubic monomials (point
values for MLS and cell averages for k-exact). The quadratic remainder about
the face centre is reproduced exactly, and the patch mean removes it. The
estimate is the least-squares solution:

```text
T = pinv(M) * [secondDerivative_k - mean_S(secondDerivative)]
deltaGradient = -scale * faceBeta * T
```

`faceBeta * pinv(M)`, with the patch mean folded in, is stored per face.
It acts directly on the cell second derivatives. Construction fails if M has
rank below ten, and the maximum condition number is reported. For cubic
fields the corrected face traction is exact with the original face rule.
Boundary faces use the same patch definition with the inward pair.

### Cell-fit correction for p=2

`curvatureCorrectionRecovery cellFit` estimates the ten third derivatives of
each cell with an unweighted least-squares cubic fit over that cell's own
reconstruction stencil (point values for MLS, cell averages for k-exact,
coordinates scaled by the stencil radius). The fit is exact for cubic
fields, so no response calibration is needed. Construction fails if the
stencil cannot support the 20-term cubic basis, and the maximum condition
number is reported; the default p=2 stencils are large enough on hex and tet
meshes. At each face the owner and neighbour estimates are averaged; a
boundary face uses the owner:

```text
T_f = (T_P + T_N)/2                 (internal)
T_f = T_P                           (boundary)
deltaGradient = -scale * faceBeta * T_f
```

Unlike `facePatch`, the residual at a face depends only on the owner and
neighbour cell stencils, the same footprint as the alpha term, rather than
on the stencils of about twelve cells.

### Q term: face-integration error with cellFit

For a cubic displacement the exact traction varies quadratically over a
face, and a rule exact only for linear functions integrates it with an
O(h^2) error that no correction evaluated at the quadrature points can
remove. For any rule exact for linear functions, that error is

```text
∫ t dA - Σ_q w_q t(x_q) = ½ H_t : (J_exact - J_rule)
J_exact = ∫ (x - c)(x - c) dA,   J_rule = Σ_q w_q (x_q - c)(x_q - c)
```

with c the face centre and H_t the (constant) Hessian of the traction, which
is linear in the third derivatives T. With `cellFit` the per-face tensor
`M = ½ (J_exact - J_rule)/A` is stored (`faceMoment_`), and the face-average
gradient increment gains `M_jk T_ijk`, which the stress law turns into
traction like the β term. J_exact comes from the fan of triangles about the
face centre; J_rule from the actual quadrature points and weights, so the
term vanishes automatically for a quadratic-exact rule.

### Alpha compensation with facePatch and cellFit

The high-order `alpha` stabilisation adds
`alpha*impKf*(uN_f - uP_f)/|n.d|` on internal faces and
`alpha*impKf*(u_b - uP_f)/|n.dL|` on fixed-value patches. Here `uP_f` and
`uN_f` are the owner and neighbour quadratic values at the face centre. For a
cubic field this jump is nonzero, so alpha adds an O(h^2) face flux and limits
the corrected scheme to second order. `alphaStab` itself is not modified.
With `facePatch` or `cellFit`, the correction also adds the negative of the
jump's cubic part:

```text
gamma_c[m] = reconstructed face-centre value of cubic monomial m from cell c
jump_f     = sum_m (gamma_N[m] - gamma_P[m]) * T[m]    (internal)
jump_f     = sum_m (-gamma_P[m]) * T[m]                (fixed-value boundary)
deltaTraction += -scale * alpha * impKf * jump_f / |n.d|
```

T is the face-patch or cell-fit estimate above. For `facePatch` the
coefficients are again stored per face and act on the cell second
derivatives; for `cellFit` the jump response is contracted with `T_f`. The
solver passes alpha only when the momentum stabilisation type is `alpha`.
Other boundary patches and `twoCell` get no compensation. For cubic fields,
alpha plus the correction is then exact on every face.

### Response correction retained for p=1

Reconstructed derivatives themselves have errors on the first omitted
polynomial degree. A geometry-only response matrix corrects this bias:

```text
raw derivatives = R * exact derivatives   (on the omitted-degree polynomials)
estimated derivatives = inverse(R) * raw derivatives
```

`R` is 6-by-6 for p=1. Its rank is checked, and the
maximum condition number is reported. It is distinct from beta.
It is not constructed or applied for p=2. Fixed-size derivative storage has
capacity for ten components; p=1 uses only six.

## One beta row per face

For p=1, use the six factorial-scaled quadratic coordinate terms:

```text
rx^2/2, rx*ry, rx*rz, ry^2/2, ry*rz, rz^2/2
```

For p=2, use the ten cubic terms:

```text
rx^3/6, rx^2*ry/2, rx^2*rz/2, rx*ry^2/2, rx*ry*rz,
rx*rz^2/2, ry^3/6, ry^2*rz/2, ry*rz^2/2, rz^3/6
```

For each original face quadrature point, `r = x - xq`. Applying the
original gradient coefficients to one of these polynomials gives its
gradient error at that point, because its exact gradient at r=0 is zero.
Boundary-data coefficients are included where the reconstruction uses them.

Beta is a vector for each coordinate term. Calculating it requires
coordinate moments and existing coefficients, not another polynomial fit.

During geometry setup, accumulate:

```text
faceBeta[m] = sum_q weight[q] * beta[q][m] / faceArea
```

Only this face-average row is stored. Runtime evaluation does not loop
over face quadrature points to apply the correction.

For moving least squares, polynomial samples are point values. For
k-exact, they are cell averages. The helper owns an auxiliary quadrature
object of cell degree p+1 to evaluate the required geometric moments.
This does not replace or raise the reconstruction's quadrature, and does
not change analytical-source integration or alpha stabilisation.

## Gradient increment and elastic traction

For each displacement component u:

```text
deltaGrad(u)[face] = -sum_m faceBeta[m] * estimatedDerivative[m]
```

For p=1, use the mean of owner and neighbour derivative estimates at internal
faces and the owner estimate at physical boundary faces. For p=2, use the
single directional estimate from the face's two-cell pair, the face-patch
estimate, or the cell-fit estimate, as described above.
Combining the three displacement components gives the tensor deltaGrad(D).

For constant-property isotropic linear elasticity:

```text
deltaSigma = 2*mu*symm(deltaGrad(D))
           + lambda*trace(deltaGrad(D))*I

deltaTraction = scale * (normal & deltaSigma)
```

The solver uses mu and lambda from the actual `linearElastic` law.
Multiplication by only `impKf` would not represent this tensor relation.

The face-average conversion is valid here because the elastic law is linear
with constant moduli and the solver uses one normal per face. It is not a
general nonlinear-material correction.

## Verification and interpretation

Build `applications/test/quadraticGradientCorrection`, then run:

```sh
Test-quadraticGradientCorrection -case /path/to/three-dimensional-case
```

The utility reads p=1 or p=2, tests both methods, and sets analytical
fixed-displacement data in memory. For monomials through degree p+1 it checks:

- The directional tensor completion agrees with an independent cubic identity
  for axis-aligned and rotated pairs, including transverse modes and pair
  reversal.
- For p=1, corrected traction equals the quadrature-weighted analytical
  traction through degree 2, using the original face quadrature rule.
- For p=2, traction remains exact through degree 2. Cubic errors before and
  after correction are reported as diagnostics, not exactness assertions.
  With `facePatch` or `cellFit` the reported corrected cubic error is at
  rounding level.
- Scales 0, 0.5 and 1 give zero, half and full additions.
- Gradients, face values and both cell/face quadrature are unchanged when
  the option is enabled.

This tests the stated polynomial properties, not elimination of the
original face-quadrature integration error. It deliberately does not expect
`fGrad()` or alpha face values to reproduce an additional polynomial degree.

Earlier MMS results for corrections embedded in gradients and alpha values
belong to a different algorithm and must not be attributed to this variant.
Nonzero beta or exact polynomial traction alone does not prove improved
global convergence or solver robustness.

The p=2 two-cell MMS experiment on six hex and six tet meshes converged for
both reconstruction methods with scale 1. It did not restore third order:
finest-three displacement slopes were approximately 1.97 (MLS) and 1.98
(k-exact) on hex, and 2.58 and 2.50 on tet. Finest hex errors exceeded the
uncorrected baseline. These observations do not isolate whether transverse
omissions, second-derivative errors, boundary treatment or other discretisation
terms dominate. Keep the correction disabled for normal use while investigating.

With `facePatch`, alpha=0 and the original face rule, k-exact on six
structured hex meshes gave pairwise displacement orders of 4.4 to 4.7. The
finest error was 6.2 times smaller than the uncorrected alpha=0.1 baseline.
MLS required thousands of linear iterations per Newton step from the second
mesh onwards and did not converge on meshes 3 to 6.

With alpha 0.1 and the compensation above, pairwise orders were 4.1 to 4.6
for k-exact and 4.1 to 4.7 for MLS. MLS converged in at most 28 linear
iterations per Newton step. Without the compensation, both fell towards
second order (last pairs 2.00 and 2.40).

`cellFit` with alpha 0.1 and compensation matches `facePatch` on hex
(finest-three slopes 4.5, errors within 5%) with fewer linear iterations.
On tets `facePatch` needs hundreds of linear iterations per Newton step and
k-exact fails at the 1000-iteration limit on five of six meshes, so
`cellFit` is the recommended p=2 recovery.

Six meshes were not enough to see the asymptotic order. On eleven meshes
(hex up to 29,791 cells, tet up to 148,567; MLS, alpha 0.1) the corrected
scheme with the default face quadrature order `p-1 = 1` falls back to
second order: finest-three slopes 1.90 (hex) and 2.51 (tet). The
one-point-per-triangle rule integrates the quadratic traction of a cubic
field with an O(h^2) error, which the point corrections cannot remove. With
`faceQuadratureOrder 2` the slopes are 4.04 (hex) and 4.02 (tet), against
4.03 and 3.94 for p=3, with finest errors 1.24 and 1.68 times those of p=3
at 75% and 96% of the p=3 run time. With the Q term and the default
one-point rule they are 4.09 and 4.05, with finest errors 1.16 and 1.21
times those of p=3 at 42% and 77% of the p=3 time (29,791 hex cells: 23 s
against 55 s; 148,567 tets: 97 s against 125 s, each run alone). Stress
stays second order in all cases. Use `curvatureCorrectionFaceIntegration
true` with `cellFit`; k-exact p=3 (23 s and 109 s, third-order stress)
remains the reference to compare against.
