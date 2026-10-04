# High-order helpers: remaining work

This file records the implementation and validation work that remains after
introducing the common least-squares reconstruction interface and the
cell-average k-exact scheme. It was audited against the current branch on
2026-09-04. Every item below still requires implementation or focused
validation.

## Recommended next sequence

1. Reject `highOrderJacobian true` in parallel until all processor-face
   Jacobian contributions are implemented.
2. Implement the k-exact mechanical and common `alphaStab` processor-face
   Jacobian paths, including complete PETSc matrix preallocation.
3. Extend `Test-alphaStabJacobian` with processor-face coverage for both
   reconstruction schemes.
4. Restore the two-processor analytical-Jacobian regressions after the focused
   Jacobian test passes.
5. Add explicit solid-model and constitutive-law capability checks so an
   unsupported high-order selection fails during setup.

## Parallel analytical Jacobian

- [ ] Assemble the `alphaStab` derivative on processor faces for both
  `movingLeastSquares` and `kExactLeastSquares`. This requires both sides of
  the face-centre value reconstruction, their global cell addressing and the
  correct equal-and-opposite row contributions.
- [ ] Extend `hofvm::initialiseJacobian()` preallocation with all local and
  remote columns used by the processor-face `alphaStab` contributions.
- [ ] Add processor-face gradient-coefficient addressing for the k-exact
  mechanical Jacobian. The explicit k-exact `fGrad()` currently exchanges
  already evaluated owner-side gradients, whereas matrix assembly needs the
  coefficients and global cell IDs from both processor sides. The MLS
  processor-face coefficient path already exists.
- [ ] Until both processor-face Jacobian paths are complete, reject
  `highOrderJacobian true` in parallel instead of silently assembling an
  incomplete matrix. In particular,
  `insertAlphaStabIntoPETScMatrix()` currently returns immediately in
  parallel.
- [ ] Extend `Test-alphaStabJacobian` with an isolated processor-face
  Jacobian-vector-product check for both reconstruction schemes.
- [ ] After the focused tests pass, add two-processor
  `highOrderJacobian-movingLeastSquares-parallel` and
  `highOrderJacobian-kExactLeastSquares-parallel` regression variants to
  `plateHole` and `sphericalCavity`.

## Stabilisation

- [ ] Add `Test-alphaStabJacobian` to a serial regression for both
  reconstruction schemes. The utility exists but is not currently called by
  a tutorial regression script.
- [ ] Add the previously planned polynomial-preservation test for
  `alphaStab`. A continuous polynomial reproduced by the selected scheme
  should give a zero reconstructed jump on internal, processor and compatible
  fixed-value faces. Test scalar and vector fields independently; the current
  Jacobian test alone cannot detect an error shared by the residual and its
  derivative.

## Boundary conditions

- [ ] Add a focused symmetry polynomial test for scalar and vector `fGrad()`.
  The `plateHole` solver regression exercises symmetry, but it does not prove
  polynomial exactness of the mirrored reconstruction.
- [ ] Validate fixed-value boundary cells whose reconstruction stencils also
  cross processor boundaries. The focused fixed-displacement polynomial test
  currently runs only in serial.
- [ ] Either support or explicitly reject non-processor coupled patches, such
  as cyclic patches, in the high-order residual. The high-order Jacobian
  already rejects them, but k-exact `fGrad()` currently only implements the
  processor coupled-patch path.
- [ ] Decide whether other vector patch fields for which `fixesValue()` is
  true need quadrature-point evaluation. Currently only
  `fixedDisplacement` and `fixedRotation` are supported, and spatially varying
  derived conditions must override `evaluateQuadrature()`.
- [ ] Implement scalar fixed-value `patchFaceQuadValues()` before using
  k-exact face gradients for pressure. Its current implementation calls
  `NotImplemented`.

## Pressure and mixed formulation

- [ ] Implement and validate the complete high-order pressure residual and
  Jacobian path before enabling a mixed displacement-pressure formulation.
  Having a pressure reconstruction object does not by itself make the mixed
  solver high order.
- [ ] Define and test pressure boundary semantics, including prescribed-value
  quadrature data and coupled displacement-pressure Jacobian blocks.
- [ ] When a pressure-solving model selects high order, require and validate
  `highOrderCoeffs.pressure` during construction rather than failing later on
  the first call to `pressureLeastSquares()`.

## Solver and constitutive-law scope

- [ ] Add explicit capability checks for solid models using
  `highOrderResidual` or `highOrderJacobian`. The implemented high-order
  residual paths are currently confined to `linGeomTotalDispSolid`,
  `nonLinGeomTotalLagTotalDispSolid` and `nonLinGeomUpdatedLagSolid`.
- [ ] Add an explicit capability check for quadrature-point stress evaluation.
  It is currently implemented only for a single material using
  `linearElastic`, `StVenantKirchhoffElastic` or `neoHookeanElastic`; other
  laws reach the base `notImplemented()` function. Record restrictions such
  as non-zero initial stress and pressure-displacement material variants.
- [ ] Document whether the nonlinear high-order Jacobians are intended as
  approximate isotropic tangent operators or extend them to use the complete
  consistent material tangent. Add a focused full-mechanical-residual
  finite-difference Jacobian check for internal, fixed-value, symmetry and
  processor-face contributions.
- [ ] Audit pointwise Taylor extrapolations that are actually reached by the
  supported high-order solid models and boundary conditions. Replace only
  those whose point-value assumption conflicts with a k-exact cell-average
  unknown.
- [ ] Document or guard the use of `movingLeastSquares` in transient cases:
  its stored unknown is a cell-centre point value, whereas the standard time
  scheme treats a `volField` value as a cell average.

## Mesh changes and transient validation

- [ ] Add a transient manufactured-solution test demonstrating that the
  k-exact cell average enters the time scheme directly without cell
  quadrature at every time step.
- [ ] Add updated-Lagrangian and moving-mesh regressions that verify cache
  invalidation, reconstruction after geometry changes and persistence of the
  required old-time quadrature kinematics.

## Geometry and focused reconstruction tests

- [ ] Resolve the cell-quadrature volume and non-zero first-central-moment
  errors observed on the general-polyhedral spherical-cavity mesh.
- [ ] Add direct scalar and vector polynomial tests for
  `faceCentreValues()` for both reconstruction schemes, including processor
  faces. These tests also provide the cleanest basis for the planned
  polynomial `alphaStab` test.
- [ ] Add direct scalar and vector polynomial tests for `valueAtPoint()` for
  both reconstruction schemes.
- [ ] Reproduce and diagnose the reported full third-order spherical-cavity
  MPI solver failure before the first SNES norm. If the current two-processor
  p=3 regression no longer reproduces it, record the successful revision and
  remove this item.
- [ ] Keep `valueAtPoint()` for sparse probe evaluation. Add a batched
  `valuesAtPoints()` interface only if a real caller needs many evaluations of
  the same field, since the current function exchanges remote stencil values
  on every call.

## Regression, compatibility and documentation

- [ ] Add a serial `highOrderJacobian` regression for
  `movingLeastSquares`; the current cantilever regression selects only
  `kExactLeastSquares` for this variant.
- [ ] Repeat the latest boundary, processor and solver tests with supported
  OpenFOAM.org versions. For foam-extend, where high-order reconstruction is
  deliberately disabled, retain build compatibility and verify the explicit
  skip/error paths instead of attempting high-order numerical tests.
- [ ] Update the higher-order README to match the current implementation. It
  still describes symmetry and the exact serial `alphaStab` Jacobian as
  future work and contains a duplicated `Creating and accessing a
  reconstruction` heading.
- [ ] Complete the per-function `kExactLeastSquares/README.md` documentation
  for the remaining lifecycle, lazy-accessor and processor-communication
  functions. The detailed derivations currently cover the principal cell and
  face coefficient functions, but not every class member.
- [ ] Record the existing structured, unstructured and tetrahedral
  manufactured-solution convergence results for polynomial orders one to
  three in a reproducible form suitable for review. Keep the unresolved
  general-polyhedral result separate.
