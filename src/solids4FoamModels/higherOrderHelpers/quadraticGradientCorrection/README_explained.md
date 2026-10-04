# The p=2 order problem and the cellFit correction, explained

This is a plain-language account of the investigation in
`CORRECTION_INVESTIGATION.md` and of the code in this directory. The
technical details are in `README.md`; this file explains the ideas.

## The problem in one paragraph

With a quadratic reconstruction (p=2) the displacement converges at second
order, not third, although a cubic reconstruction (p=3) gives fourth order
and a linear one (p=1) gives second. The reason is that p=2 is an even
order. Its leading error is an odd term (the cubic one), and odd terms do
not cancel between the two sides of a face. The only cure is to give the
face flux the cubic information it lacks. That is what the correction does,
and `cellFit` is the cheapest way we found to do it.

## Why an even order loses one order

Think of one face between cell P (owner) and cell N (neighbour). Each cell
reconstructs the displacement with a polynomial and evaluates it, or its
gradient, at the face.

- A degree-p polynomial fitted to the data reproduces every term of the true
  solution up to degree p. The first thing it gets wrong is the degree p+1
  term.
- For p=1 the missed term is quadratic, an even function of the distance
  from the face. The owner and the neighbour make the same mistake, so the
  mistake cancels in the shared flux. The flux is better than expected and
  the scheme gains an order: second order from a linear reconstruction.
- For p=3 the missed term is quartic, also even. Same cancellation, so p=3
  gives fourth order.
- For p=2 the missed term is cubic, an odd function. The owner's mistake and
  the neighbour's mistake have opposite signs, so they add up instead of
  cancelling. The flux error is O(h^2), and a finite-volume solution has the
  order of its flux error: second order.

In short: p=2 has a systematic O(h^2) flux error proportional to the third
derivatives of the displacement, and it survives.

The alpha stabilisation has the same weakness. It adds a traction
proportional to the jump between the owner's and neighbour's face values.
For a cubic field these two quadratic face values differ, so alpha adds an
O(h^2) traction of its own. In the 1-D script, uniform mesh, the alpha jump
for p=2 is 1.79 h^3 at every mesh size, while for p=3 it is zero to rounding.
Alpha is harmless whenever it is not the largest error, which is why it
never showed up at p=1 or p=3 or in uncorrected p=2.

## How the correction works

The reconstruction is linear in the data, so its error on a cubic field can
be computed once from the geometry. For each face we apply the actual
gradient stencil to the ten cubic monomials (x^3/6, x^2 y/2, ..., z^3/6).
The result, called beta, is a 3-by-10 table per face: how much gradient
error each unit third derivative produces. Beta depends only on the mesh and
the stencil weights, never on the solution.

At run time the correction needs an estimate of the ten third derivatives T
of the current displacement at the face. Then

```text
gradient error at the face  =  beta * T
correction traction         = -(stress law applied to beta * T)
```

is added once per face as an extra traction. With `cellFit` the same T also
removes the cubic part of the alpha jump, using a similar geometry-only
table. Nothing else changes: `fGrad`, the face values, the alpha model and
the quadrature are untouched. Setting the scale to zero recovers the
original residual exactly.

The whole difficulty is estimating T well enough. It only needs O(h)
accuracy, because it multiplies an O(h^2) coefficient, but it must be a
consistent estimate of all ten components.

## The three ways of estimating T

- `twoCell` (the original attempt): difference the owner's and neighbour's
  second derivatives along the line between them. In 1-D this is the whole
  answer. In 3-D it only sees the six third derivatives with at least one
  index along that line; the four purely transverse ones are invisible. The
  audit showed that no combination of the pair's data can cancel the cubic
  error (rank 6 of 10), so this variant could never work in 3-D.
- `facePatch`: use the second derivatives of the owner, the neighbour and
  all their face neighbours (about twelve cells) and solve a small
  least-squares problem for all ten T. This is exact for cubic fields and
  gave the first third-order results. Its drawback is cost: each of those
  twelve cells has its own 65-cell stencil, so the residual at a face
  depends on a "stencil of stencils". On tets the linear solver needed
  hundreds of iterations per Newton step and often failed.
- `cellFit` (recommended): fit a cubic polynomial to each cell's own p=2
  stencil, read off the ten third derivatives, and average owner and
  neighbour at the face. The next section explains this step by step.

## cellFit step by step

What is needed at each face is T, the ten third derivatives of the
displacement: u_xxx, u_xxy, u_xxz, u_xyy, u_xyz, u_xzz, u_yyy, u_yyz,
u_yzz, u_zzz. Beta is known from the geometry, so T is the only unknown.

Step 1: fit a cubic to each cell's own stencil. Take one cell P. The p=2
reconstruction already has a stencil for it: the cell itself plus about 65
neighbours, with the displacement values in those cells (point values for
MLS, cell averages for k-exact). `cellFit` fits a cubic polynomial through
those same values by least squares:

```text
u(x) ~ a0 + a1 x + a2 y + a3 z
     + a4 x^2/2 + ... + a9 z^2/2                 (6 quadratic terms)
     + a10 x^3/6 + a11 x^2 y/2 + ... + a19 z^3/6 (10 cubic terms)
```

That is 20 unknowns fitted to about 65 data points. The coordinates are
measured from the cell centre, and each monomial is divided by its
factorial on purpose: with that scaling every coefficient is a derivative
of the fitted polynomial at the cell centre.

Step 2: read off the ten third derivatives. This means nothing more than
taking the last ten coefficients: a10 is u_xxx, a11 is u_xxy, and so on up
to a19, which is u_zzz. These ten numbers are T_P, the estimate for cell P.

Two things make this cheap:

- The fit is linear in the data, so the least-squares solution is a fixed
  set of weights per cell, `T_P[m] = sum_j w[P][j][m] u_j` over the stencil
  cells. The weights are computed once from the geometry at start-up
  (`makeCellFitCoeffs`). At run time each cell does a 65-by-10 weighted sum.
- The fit is exact for cubic fields: if the true displacement is a cubic
  polynomial, the fit reproduces it and T is exact. For a general smooth
  field the error in T is O(h), which is all the correction needs, because
  T multiplies an O(h^2) coefficient.

Step 3: get T at the face. The face between P and N has two estimates, T_P
and T_N. `cellFit` takes their average, `T_f = (T_P + T_N)/2`. A boundary
face has only one cell, so `T_f = T_P`. The average is the symmetric
choice; either estimate alone would also be O(h) accurate.

Step 4: apply the correction. The extra traction is the stress law applied
to `-beta * T_f`, plus the alpha compensation built from the same T_f.

Why the footprint matters. The footprint of a face is the set of cells its
flux depends on. In the uncorrected p=2 scheme the alpha term already makes
the flux depend on the owner's and the neighbour's reconstruction stencils.
`cellFit` uses exactly those two stencils, so the corrected flux depends on
no extra cell. `facePatch`, by contrast, used the curvatures of about
twelve cells, each with its own 65-cell stencil, so a flux depended on
several hundred cells. That widening is what made it hard to solve on tets.

Why the stencils are large enough. A cubic fit has 20 unknowns, so it needs
at least 20 well-spread data points, and comfortably more for good
conditioning. A minimal p=2 stencil (10 unknowns) might have only 20 to 30
cells. The stencils here have about 65 because about 55 extra cells are
added for conditioning. The check in
`cubic-consistency-audit.I1sEzZ/stencilCheck` ran the cubic fit on every
cell of every hex and tet mesh: full rank everywhere, condition numbers
between 20 and 300.

In one sentence: `cellFit` is a p=3-style fit done on the p=2 stencil, used
only to supply the ten third derivatives to the correction, while the
reconstruction, the fluxes and alpha stay p=2.

`cellFit` is selected in the displacement reconstruction dictionary:

```text
curvatureCorrection         true;
curvatureCorrectionScale    1;
curvatureCorrectionRecovery cellFit;
faceQuadratureOrder         2;
```

## What the results show

Manufactured solution, linear elastic, alpha 0.1, MLS reconstruction,
eleven meshes per family (hex up to 29,791 cells, tet up to 148,567).
Fitted displacement order on the three finest meshes, with the finest-mesh
error in metres and the solver time in seconds:

| Variant | Hex | Tet |
| --- | ---: | ---: |
| p=2 | 2.26 / 1.4e-10 / 12 | 2.15 / 6.2e-11 / 49 |
| p=2 + cellFit, face order 1 | 1.90 / 1.4e-11 / 23 | 2.51 / 1.8e-11 / 96 |
| p=2 + cellFit, face order 2 | 4.04 / 3.2e-12 / 41 | 4.02 / 6.5e-13 / 119 |
| p=3 | 4.03 / 2.6e-12 / 56 | 3.94 / 3.8e-13 / 137 |

- With the default face quadrature (one point per triangle) the correction
  cuts the error a lot but the order returns to two on fine meshes. The
  corrected flux is exact for cubics at the quadrature points, but the
  exact traction of a cubic field varies quadratically over the face, and
  the one-point rule integrates that with an O(h^2) error.
- With a quadratic-exact face rule (`faceQuadratureOrder 2`) the corrected
  p=2 is fourth order on both mesh families, like p=3, with errors 1.2 to
  1.7 times those of p=3 at 74% to 87% of the p=3 time. It is fourth order
  rather than third because the remaining error after the cubic term is the
  quartic one, which is even and cancels across the face.
- Stress stays second order because only the face traction is corrected.
- Six meshes were not enough to see any of this: on the coarse meshes the
  corrected scheme showed 4.5 on hex and 2.9 on tets regardless of the face
  rule.

## What was ruled out along the way

Each of these was tested and did not explain the order loss:

- Boundary cells dominating the error norm (interior-only norms gave the
  same order).
- A bug in beta or in the traction assembly (independent recomputation
  agrees to 1e-10).

What did explain it: the two-cell estimate cannot see four of the ten third
derivatives; alpha adds its own O(h^2) error once the reconstruction's
error is removed; and the one-point face rule adds a third O(h^2) error
that only shows once the other two are gone and the mesh is fine enough.
The face rule was first dismissed because raising it changed nothing on the
six coarse meshes.

## What is not done

- Only serial, constant-property `linearElastic`, `solvePressure false`
  and a matrix-free Jacobian are supported.
- Stress is not corrected and stays second order; that was the intent.
- Six coarse meshes cannot fix an asymptotic order; tet pairwise orders
  scatter between 2.3 and 3.6.
- No timing on large meshes has been done. Whether `cellFit` p=2 or plain
  p=3 is cheaper at a given accuracy is still open.
