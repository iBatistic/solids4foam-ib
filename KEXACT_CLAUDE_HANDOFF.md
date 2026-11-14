# k-Exact High-Order Framework: Claude Handoff

## Purpose

This document summarises the design decisions, implementation history, current
code structure, validation utilities, and unfinished work for the solids4foam
high-order least-squares framework. It is intended to be supplied to Claude at
the beginning of a new session.

Claude should treat the current source code as authoritative. This document
explains why the code has its present structure and identifies areas that are
known to remain incomplete.

## Current repository state

At the time this handoff was written:

- repository: `/Volumes/OpenFOAM/work/solids4foam-ib/solids4foam-ib`;
- branch: `feature-highOrderUpdates`;
- `HEAD`: `b200718f`;
- `feature-highOrderUpdates`, local `development`, and `origin/development`
  point to the same commit;
- the worktree was clean before this handoff file was added;
- the original k-exact feature was merged as commit `202c50d7`,
  `Feature: k-exact least squares (#398)`.

The original implementation conversation used an older checkout named
`solids4foam-ib-2412` and branch `feature-kExactLeastSquares`. Do not assume
that old paths or branch state still apply. The implementation has subsequently
received foam-extend support, solid-model integration, stencil changes,
quadrature corrections, and explicit wedge-mesh rejection.

## Original motivation

The original `movingLeastSquares` implementation interprets a stored
cell-centre field value as the point value at `mesh.C()[cellI]`. It reconstructs
derivatives by building weighted least-squares systems around the requested
point, including face quadrature points.

That interpretation is inconvenient for transient problems. A point-valued
unknown must be integrated over every cell to form the time term. Repeating
that volume integration at every time step is expensive.

The k-exact formulation instead stores the cell average

$$
\overline{u}_P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}u(\boldsymbol{x})\,\mathrm{d}\Omega.
$$

The time discretisation can therefore use the stored unknown directly. Spatial
derivatives and values away from the centroid are obtained from a mean-free
polynomial reconstruction that preserves the owner-cell average.

## Naming and architecture decisions

The agreed names are:

- abstract runtime-selectable interface: `leastSquaresScheme`;
- abstract per-cell-polynomial implementation base: `cellLeastSquares`
  (not run-time selectable);
- moving evaluation-point reconstruction: `movingLeastSquares`;
- point-value, one-polynomial-per-cell reconstruction: `leastSquares`;
- cell-average, moment-corrected reconstruction: `kExactLeastSquares`;
- shared stencil class: `leastSquaresStencil`;
- mesh-object owner: `leastSquaresReconstruction`;
- solid-model accessors: `displacementLeastSquares()` and
  `pressureLeastSquares()`;
- face-gradient coefficient builder: `calcFaceGradCoeffs()`.

`pointwiseLeastSquares`, `movingLeastSquaresStencil`, `displacementMLS`, and
`calcFaceCoeffs` are obsolete names and should not be reintroduced.

The principal ownership hierarchy is:

```text
fvMesh
└── leastSquaresReconstruction              mesh object, lazy owner
    ├── displacement leastSquaresScheme     selected from dictionary
    └── pressure leastSquaresScheme         selected from dictionary
        ├── movingLeastSquares              point-valued unknowns
        └── cellLeastSquares                one polynomial per cell
            ├── leastSquares                point-valued unknowns
            └── kExactLeastSquares          cell-average unknowns
```

`leastSquaresScheme::cellAverageUnknowns()` tells whether a scheme stores
cell averages (`kExactLeastSquares`) or point values (`movingLeastSquares`,
`leastSquares`). Code outside the reconstruction classes must use it, never
the type name or `isA<kExactLeastSquares>`.

`leastSquaresReconstruction` owns one lazy reconstruction per field role. This
allows displacement and pressure to select different schemes while avoiding a
high-order object stored directly in `solidModel`.

Typical access is:

```cpp
const leastSquaresReconstruction& reconstructions =
    leastSquaresReconstruction::New(mesh);

const leastSquaresScheme& displacement =
    reconstructions.scheme
    (
        "displacement",
        includePatchInStencils,
        displacementDictionary
    );
```

Solid solvers normally use `solidModel::displacementLeastSquares()` or
`solidModel::pressureLeastSquares()` rather than constructing the mesh object
directly.

## Important source files

- `src/solids4FoamModels/higherOrderHelpers/leastSquaresScheme/`
  contains the abstract interface and shared face-centre value evaluation.
- `src/solids4FoamModels/higherOrderHelpers/leastSquaresReconstruction/`
  contains the mesh-object owner.
- `src/solids4FoamModels/higherOrderHelpers/movingLeastSquares/`
  contains the original point-valued method.
- `src/solids4FoamModels/higherOrderHelpers/cellLeastSquares/`
  contains the shared one-polynomial-per-cell implementation (coefficient
  builds, point evaluation, face gradients, boundaries, communication).
- `src/solids4FoamModels/higherOrderHelpers/kExactLeastSquares/`
  contains the thin cell-average scheme and the detailed mathematical README
  of the shared implementation.
- `src/solids4FoamModels/higherOrderHelpers/leastSquares/`
  contains the thin point-value scheme and a README on its reduction from
  the k-exact formulation.
- `src/solids4FoamModels/higherOrderHelpers/leastSquaresStencil/`
  contains local and off-processor stencil construction and field exchange.
- `src/solids4FoamModels/higherOrderHelpers/fvMeshQuadrature/`
  contains cell and face quadrature, cell moments, and processor-face
  quadrature synchronisation.
- `src/solids4FoamModels/higherOrderHelpers/hofv/`
  contains high-order explicit and PETSc matrix operators.
- `src/solids4FoamModels/numerics/stabilisationModels/alphaStab/`
  contains the alpha stabilisation residual.

Before changing the k-exact mathematics, read
`src/solids4FoamModels/higherOrderHelpers/kExactLeastSquares/README.md`. It
documents the equations and the exact ordering of coefficient arrays.

## k-Exact formulation implemented in the code

For owner cell $P$, coordinates are measured from the volume centroid:

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P.
$$

The mean-free reconstruction through cubic order is

$$
\begin{aligned}
u_P^R(\boldsymbol{x})={}&\overline{u}_P
+\boldsymbol{r}\boldsymbol{\cdot}(\nabla u)_P \\
&+\frac{1}{2}(\nabla\nabla u)_P:
\left(
\boldsymbol{r}\otimes\boldsymbol{r}-\boldsymbol{M}^{(2)}_P
\right) \\
&+\frac{1}{6}(\nabla\nabla\nabla u)_P:
\left(
\boldsymbol{r}\otimes\boldsymbol{r}\otimes\boldsymbol{r}
-\boldsymbol{M}^{(3)}_P
\right).
\end{aligned}
$$

Subtracting the owner-cell moments makes every non-constant basis term have
zero mean over $P$. Consequently, the average of the reconstructed polynomial
is exactly $\overline{u}_P$.

The implementation supports polynomial orders 1, 2, and 3. In 2-D, all terms
containing $z$ are omitted, and the empty direction must be $z$. Three-
dimensional meshes are supported. Wedge/axisymmetric meshes are deliberately
rejected.

`calcCellCoeffs()` builds a weighted least-squares system using differences of
cell averages,

$$
\overline{u}_N-\overline{u}_P,
$$

and volume-averaged mean-free monomials. It solves the small dense systems with
Eigen QR and stores geometry-only coefficients for first, second, and third
derivatives. The same coefficients can then be applied to scalar fields or to
each component of a vector field.

`calcCellCoeffs()` does not incorporate boundary-condition observations. The
accepted final design keeps cell-centred derivative coefficients independent
of boundary field values. Face-specific boundary reconstruction is handled in
`calcFaceGradCoeffs()`.

## Cell moments and quadrature

First-, second-, and third-order volume-normalised central moments are owned
and lazily evaluated by `fvMeshQuadrature`, not by `kExactLeastSquares`.

The first moments should be zero when `mesh.C()` is the true volume centroid.
Second moments are required for $p \ge 2$. Third moments are only requested by
the k-exact method for $p=3$.

For a 2-D mesh, moments are normalised by cell volume as represented in the
OpenFOAM mesh, while all components involving the empty $z$ direction are
omitted or zeroed. Do not introduce a separate area-normalisation convention
without auditing OpenFOAM's extruded-volume representation.

Processor-face quadrature points and weights are synchronised in
`fvMeshQuadrature`. The lower processor rank supplies the ordering used on
both sides. Consequently, reconstruction code should not repeat point-location
matching inside `fGrad()`.

## Stored coefficient conventions

The documentation uses the following geometry-dependent coefficient notation:

- $C_x$ for first-derivative coefficients;
- $C_{xx}$ for second-derivative coefficients;
- $C_{xxx}$ for third-derivative coefficients;
- $C$ for value coefficients at an arbitrary point.

For a stencil cell $N$, derivative rows multiply the cell-average difference:

$$
(\nabla u)_P
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{x,N}^P
\left(\overline{u}_N-\overline{u}_P\right).
$$

`cellValueCoeffsAtPoint(cellI, x, coeffs)` collapses the first-, second-, and
third-derivative effects into scalar coefficients for evaluating the value at
`x`. Its output is ordered as all stencil cells followed by the owner cell.
The absolute coefficients sum to one.

`cellGradCoeffsAtPoint(cellI, x, coeffs)` collapses the same derivative rows
into vector coefficients for the gradient at `x`. It returns one entry per
stencil cell and does not append the owner. The implied absolute owner
coefficient is minus the sum of the stencil coefficients.

`valueAtPoint(field, cellID, x)` is implemented for both
`movingLeastSquares` and `kExactLeastSquares`. It evaluates only the requested
cell row once mesh-wide geometry tables have been lazily constructed. In
parallel it is a collective pattern: every rank participates in remote-field
exchange, the owning rank passes its local `cellID`, and other ranks pass `-1`.
`solidPointDisplacement` uses this pattern.

This interface is suitable for a small number of requested points. It is not a
batched API optimised for repeatedly evaluating every cell in a large loop.

## Face-gradient reconstruction

For an internal face, the face stencil merges the owner cell and stencil with
the neighbour cell and stencil. Duplicate global cell IDs are removed. At each
face quadrature point, the final gradient is the average of the two cell
reconstructions.

Cells that occur on only one side simply receive the coefficient contributed
by that side. Cells present in both reconstructions receive the sum of both
side contributions. Owner and neighbour effects are therefore fully collapsed
into one coefficient row for an internal face.

For processor faces, `kExactLeastSquares::fGrad()` currently evaluates the
owner-side reconstruction on each rank, exchanges the reconstructed gradient,
and averages the two values. It does not construct a merged processor-face
sparse coefficient row.

## Boundary treatment

### Prescribed displacement and rotation

An early exact-constraint formulation was removed because boundary constraint
matrices could be rank deficient. A representative failure was six constraints
with rank five for nine polynomial derivative unknowns.

The accepted formulation adds prescribed values at each boundary face
quadrature point as weighted least-squares observations in the face-specific
system. They are not exact algebraic constraints. The weights come from the
same selected weight function as the cell observations.

`patchFaceQuadValues()` obtains quadrature-point values from:

- `fixedDisplacementFvPatchVectorField::evaluateQuadrature()`; and
- `fixedRotationFvPatchVectorField::evaluateQuadrature()`.

This detail was essential for recovering the expected high-order convergence.
Using only a face-centre boundary value caused the cubic reconstruction to
converge incorrectly.

Known boundary data have separate addressing and coefficients from global cell
unknowns. When differentiating the residual for the Jacobian, prescribed face
values have zero derivative. Their value still appears in the explicit
residual, so the corresponding nonhomogeneous contribution is naturally on
the residual/right-hand-side side of the nonlinear equation.

### Traction and other Neumann boundaries

Physical traction/fixed-gradient boundaries use a real one-sided owner-cell
gradient reconstruction. They are not left at zero and are not treated as
prescribed displacement observations.

### Symmetry

Symmetry uses a reflected stencil, following the idea previously implemented
in `movingLeastSquares`. The face reconstruction contains physical stencil
values and their reflected counterparts. Vector values are transformed with
the Householder reflection tensor.

The `plateHole` tutorial was used during development to confirm that the
symmetry treatment produced a converged solution.

## Alpha stabilisation and analytical Jacobian

The shared `leastSquaresScheme::faceCentreValues()` interface evaluates
owner- and neighbour-side reconstructed values at face centres for both
`movingLeastSquares` and `kExactLeastSquares`. The alpha stabilisation residual
uses this common path.

`hofvm::insertAlphaStabIntoPETScMatrix()` inserts the analytical derivative of
the alpha-stabilisation term. In serial it currently covers:

- both owner- and neighbour-side reconstructions on internal faces; and
- the owner reconstruction on fixed-value boundaries.

For an internal face, the same flux contribution is inserted into the owner
row and with the opposite sign into the neighbour row. Both reconstructed face
values depend on cell unknowns, which is why the helper is called once for the
owner reconstruction and once for the neighbour reconstruction.

On a fixed-value boundary, the prescribed face value is independent of cell
unknowns. Only the derivative of the extrapolated owner reconstruction is
inserted into the matrix.

## Current configuration pattern

A displacement reconstruction is selected under the solid model's
`highOrderCoeffs` dictionary:

```text
highOrderCoeffs
{
    highOrderJacobian false;
    highOrderResidual true;

    displacement
    {
        type kExactLeastSquares;
        polynomialOrder 3;
        cellStencilExtraCells 15;

        weightFunctionCoeffs
        {
            type radiallySymmetricExponential;
            k 6;
        }

        calcConditionNumber false;
    }
}
```

For backward compatibility, `kExactLeastSquares` reads
`faceStencilExtraCells` when `cellStencilExtraCells` is absent. New cases
should use `cellStencilExtraCells`.

High-order residual and Jacobian paths currently require the `PETScSNES`
solution algorithm.

Analytical-solution function objects check the displacement reconstruction
type. They integrate the analytical displacement over each cell and divide by
`mesh.V()` for `kExactLeastSquares`; moving least squares and conventional
methods compare against the analytical point value at the cell centre.

## Validation utilities

### `Test-highOrderGrad`

Location: `applications/test/highOrderGrad/`.

The utility reads the reconstruction settings from `solidProperties` and
tests:

- scalar and vector fields;
- linear, quadratic, and cubic polynomials;
- cell-centred `grad()`;
- internal-face quadrature-point `fGrad()`;
- processor-face quadrature-point `fGrad()` in parallel;
- traction/Neumann boundary gradients;
- fixed-displacement boundary gradients in the supported serial focused test;
- both 2-D and 3-D meshes.

For k-exact, it initialises exact cell averages. For moving least squares, it
initialises point values at cell centres. Only polynomial degrees up to the
configured reconstruction order determine pass/fail.

### `Test-cellMoments`

Location: `applications/test/cellMoments/`.

The utility prints first-, second-, and third-order central moments for every
cell. It also prints:

- average and maximum first-moment magnitude;
- the sum of cell quadrature weights;
- relative error against `mesh.V()` for every cell; and
- average and maximum relative volume error.

Use this utility first when k-exact polynomial reproduction is poor on a
polyhedral mesh.

### `Test-alphaStabJacobian`

Location: `applications/test/alphaStabJacobian/`.

The utility compares the analytical matrix-vector product with a centred
finite-difference Jacobian-vector product of the explicit alpha-stabilisation
residual. Internal-face and fixed-value-boundary contributions are isolated
with face-diffusivity masks. It applies to both reconstruction types.

It currently skips when PETSc is unavailable and deliberately skips parallel
execution.

### Tutorial regression coverage

The linear-elasticity regressions include moving- and k-exact variants in
cases such as:

- `cantilever2d`;
- `cooksMembrane`;
- `plateHole`;
- `pressurisedCylinder`; and
- `sphericalCavity`.

`plateHole` and `sphericalCavity` include two-processor high-order residual
variants for both reconstruction types. `cooksMembrane` runs the serial
gradient utility, and `sphericalCavity` runs its parallel form when available.

Historical development runs showed close displacement results between moving
and k-exact reconstructions after the boundary quadrature fix. These results
were not rerun while creating this handoff, so current validation claims should
come from fresh logs rather than this historical observation.

## Parallel scaling regression after #398 (resolved 2026-11-03)

After PR #398 the parallel wall-clock time of high-order runs was flat with
the number of ranks (a friend measured about 99 s for 2, 3 and 4 ranks against
109 s in serial on a 91k-cell trigonometric manufactured-solution case). The
cause was not the face-centre evaluation of the alpha stabilisation, which
scales. It was `fvMeshQuadrature::synchroniseProcessorFaceQuadrature()`,
added in #398: the lower rank posted buffered sends only and then went on to
`calcFaceCoeffs()`. A buffered MPI send is only progressed by the sender's
next MPI call, so the higher rank sat in its blocking receive until the lower
rank finished its coefficient build. This serialised the one-time face
coefficient build rank by rank. The function now packs the quadrature data
per neighbouring processor and exchanges it collectively: `Pstream::exchange`
on foam-extend, which has no `PstreamBuffers`, and `PstreamBuffers` with
`finishedSends()` before any receive on OpenFOAM.com and OpenFOAM.org. On
foam-extend each `Pstream::exchange` call takes its own `labelListList` size
table: foam-extend fills the table with a `listEq` combine, so a table reused
from a previous call keeps the stale entries of the other ranks, the receive
is posted with the wrong byte count and `MPI_Waitall` fails with
`MPI_ERR_TRUNCATE` (CI run on PR #662, `plateHole` and `sphericalCavity`
parallel high-order variants). The file compiles against OpenFOAM v2412,
OpenFOAM 9 and foam-extend 4.1; the `plateHole` high-order variants were run
on v2412 and foam-extend 4.1.

Measured on the Apple M4 Pro with OpenFOAM v2412, 91 125 hexahedra,
`movingLeastSquares` p=3, 15 extra stencil cells, hypre preconditioner:

| ranks | before fix | after fix |
|-------|------------|-----------|
| 1     | 24-28 s (30k cells) | 80.3 s |
| 2     | 20.5 s (30k cells)  | 42.6 s |
| 4     | 18.4 s (30k cells)  | 22.2 s |
| 8     | n/a                 | 15.4 s |

Before-fix numbers are for the 30k-cell mesh, after-fix for the 91k-cell
mesh; the error norms are identical for all rank counts. Do not reintroduce
per-patch blocking send/receive pairs in which one side only sends.

`kExactLeastSquares` on the 30k-cell mesh with 60 extra stencil cells (the
value the sphericalCavity tutorial uses for it) scales the same way: 23.7,
13.4 and 7.9 s at 1, 2 and 4 ranks with identical error norms. With 15 extra
cells the k-exact linear solve diverges already in serial, which is a case
setting and unrelated to the fix. The sphericalCavity regression script
passes for all six approaches; the plateHole script needs bash 4 (associative
arrays), so its four high-order variants were run by hand at 2 ranks and are
within its tolerances.

The benchmark case's `manufacturedSolution` boundary condition does not write
its `E`, `nu`, `ax`, `ay`, `az` entries, so `decomposePar` produces fields
that fail to read; re-insert them in `processor*/0/D` or fix its `write()`.

## Confirmed unfinished work

### Parallel analytical high-order Jacobian

This is the main unfinished implementation task.

`hofvm::insertAlphaStabIntoPETScMatrix()` immediately returns when
`Pstream::parRun()` is true. The processor-face alpha-stabilisation derivative
is therefore not assembled.

There is a second k-exact-specific issue that must be considered at the same
time: k-exact processor-face `fGrad()` exchanges already reconstructed
gradients, while `makeFaceGradStencil()` deliberately leaves coupled-face rows
empty. The PETSc Jacobian needs explicit sparse global-cell addressing and
coefficients. Removing only the early return in the alpha helper is not enough
to establish a correct parallel k-exact analytical Jacobian.

Before enabling `highOrderJacobian true` in parallel regressions, complete and
verify all of the following:

1. Define sparse processor-face derivative rows for both sides of a k-exact
   face, using global cell IDs.
2. Exchange the remote reconstruction addressing and coefficients, or provide
   an equivalent deterministic construction on each rank.
3. Assemble each rank's local residual row with the correct owner-side sign.
   The neighbouring rank assembles its own row.
4. Add all new local and off-processor columns to
   `hofvm::initialiseJacobian()` preallocation. The matrix has
   `MAT_NEW_NONZERO_ALLOCATION_ERR` enabled, so missing columns will fail.
5. Add processor-face alpha-stabilisation contributions.
6. Extend `Test-alphaStabJacobian` with a processor-face diffusivity mask and a
   two-rank finite-difference comparison.
7. Only after the focused test passes, add two-rank
   `highOrderJacobian` regressions for moving and k-exact reconstruction.

Until this work is complete, use `highOrderJacobian false` for parallel runs.
The usual second-order Jacobian/preconditioner can still be used with the
high-order residual.

### Additional limitations and follow-up candidates

- Cyclic coupled patches are not implemented for the high-order analytical
  Jacobian.
- The focused fixed-displacement boundary-gradient test currently runs only in
  serial.
- `valueAtPoint()` is collective in parallel and is not designed as a batched
  whole-field evaluator.
- Polyhedral meshes with inaccurate centroids or cell quadrature should be
  diagnosed with `Test-cellMoments`; do not compensate by changing only the
  moment normalisation.
- The method supports polynomial orders only from 1 through 3.
- Wedge meshes are intentionally unsupported.

## Changes that should not be repeated or reverted

- Do not move moments back into `kExactLeastSquares`; they belong to
  `fvMeshQuadrature` and are lazily cached there.
- Do not restore exact fixed-value constraints. Weighted quadrature-point
  observations are the accepted formulation.
- Do not add boundary observations to `calcCellCoeffs()` without a new design
  review. Boundary systems are face specific and currently belong to
  `calcFaceGradCoeffs()`.
- Do not treat stored k-exact unknowns as point values when computing error
  norms or time terms.
- Do not reintroduce field-wide `cellDerivatives()` merely for
  `valueAtPoint()`; the existing interface evaluates the selected cell row.
- Do not duplicate processor quadrature-point matching in `fGrad()`;
  synchronisation belongs to `fvMeshQuadrature`.
- Do not enable parallel high-order-Jacobian regression merely because a case
  reaches `End`; first establish the focused Jacobian-vector-product check.

## Literature and user derivations

The implementation was developed using the following references:

1. C. B. Sejekan and C. F. Ollivier-Gooch, "Improving finite-volume
   diffusive fluxes through better reconstruction," *Computers & Fluids*,
   139, 216-232, 2016.
2. A. Jalali and C. F. Ollivier-Gooch, "Higher-Order Finite Volume Solution
   Reconstruction on Highly Anisotropic Meshes," 21st AIAA Computational
   Fluid Dynamics Conference, 2013.
3. Ricardo Daniel Pereira da Costa, *Advanced computational methods towards
   high-performance polymer processing simulations*, PhD thesis, 2020.
4. C. Ollivier-Gooch and M. Van Altena, "A High-Order-Accurate Unstructured
   Mesh Finite-Volume Scheme for the Advection-Diffusion Equation,"
   *Journal of Computational Physics*, 181, 729-752, 2002.

User-provided local documents used during the derivation included:

```text
/Users/ivanbatistic/Downloads/1-s2.0-S0021999126004092-main (1).pdf
/Users/ivanbatistic/Downloads/main.tex
/Users/ivanbatistic/Downloads/main (1).tex
/Users/ivanbatistic/Downloads/eq_derivation.tex
/Users/ivanbatistic/googleDrive/Literature/High-order/
  A High-Order-Accurate Unstructured Mesh Finite-Volume Scheme for the
  Advection-Diffusion Equation.pdf
/Users/ivanbatistic/googleDrive/Literature/High-order/High-order-CFD/
  Costa-Nobrega/costaThesis.pdf
```

These paths may need shell quoting because they contain spaces and
parentheses.

## Recommended start for the next Claude session

Ask Claude to:

1. read this handoff;
2. inspect the current Git status and branch rather than assuming the old
   feature branch;
3. read the k-exact README and the current implementations of
   `hofvm::initialiseJacobian()`, `hofvm::laplacianIntoPETScMatrix()`, and
   `hofvm::insertAlphaStabIntoPETScMatrix()`;
4. inspect the processor-face paths in `kExactLeastSquares::fGrad()` and
   `leastSquaresStencil`;
5. propose the processor-face sparse addressing, communication, and PETSc
   preallocation design before editing code; and
6. preserve OpenFOAM.com, OpenFOAM.org, and foam-extend compatibility.

No fresh build or numerical test was run while preparing this handoff. Its
current-state statements are based on source and Git inspection.
