# Plan: `leastSquares` reconstruction scheme and `cellLeastSquares` base

Status: implemented on `feature-highOrderUpdates_02` (2026-11-10) on top of
upstream `development` 20643798, which includes #662 and #663; under
review. Revision 2 adopts the agreed hierarchy and naming:

```text
leastSquaresScheme                 common runtime interface
├── movingLeastSquares             moving evaluation-point reconstruction
└── cellLeastSquares               abstract per-cell-polynomial implementation
    ├── leastSquares               point-value, one polynomial per cell
    └── kExactLeastSquares         cell-average, moment-corrected
```

The plural spelling `LeastSquares` is used throughout. `cellLeastSquares` is
an implementation base and is not registered in the runtime-selection table.
`movingLeastSquares` stays derived directly from `leastSquaresScheme`; the
asymmetric depth is intentional.

## 1. Goal

Add a third run-time selectable scheme, `type leastSquares;`, next to
`movingLeastSquares` and `kExactLeastSquares`:

- the stored cell unknown $u_P$ is the field value at the cell centre;
- one polynomial per cell, centred at the owner-cell centre, with the constant
  fixed by $u_P$, so $u_P^R(\boldsymbol{x}_P)=u_P$;
- neighbour observations are point values at their cell centres;
- no cell moments;
- the cell polynomial is reused at cell centres, face centres, face
  quadrature points and arbitrary points.

It takes the same dictionary entries as `kExactLeastSquares`:
`polynomialOrder`, `cellStencilExtraCells` (with `faceStencilExtraCells`
accepted as today), `haloDepthScale`, `weightFunctionCoeffs`,
`calcConditionNumber`.

The class `Description` and the README must state explicitly that
`leastSquares` stores point values at cell centres, because OpenFOAM also has
a gradient scheme called `leastSquares` in a different selection table.

## 2. Mathematics: `leastSquares` is k-exact with zero positive-order moments

For owner cell $P$ and $\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P$ the
k-exact reconstruction is

$$
u_P^R(\boldsymbol{x})=\overline{u}_P
+\sum_{1\le|\alpha|\le p}\frac{\partial^{\alpha}u_P}{\alpha!}
\left(\boldsymbol{r}^{\alpha}-M_P^{(\alpha)}\right),
$$

with $M_P^{(\alpha)}$ the volume-normalised central moments of $P$. The
stencil observations are differences of cell averages, each average being
`averageMonomial()`, a binomial sum of the neighbour's central moments times
powers of $\boldsymbol{d}_{PN}$.

Setting every positive-order moment to zero,

$$
M^{(\alpha)}\;\to\;\delta_{|\alpha|,0},
$$

gives the pointwise Taylor fit:

- `averageMonomial(N, d, alpha)` collapses to the single term
  $\boldsymbol{d}_{PN}^{\alpha}$, so the cell system is
  $u_N-u_P=\sum_\alpha \boldsymbol{d}_{PN}^{\alpha}/\alpha!\;\partial^\alpha u_P$;
- the owner-cell averages (`fillOwnerAverages`) are zero for every
  non-constant term;
- `reflectedCentralMoment` already returns 1 for order 0 and 0 for order 1,
  and its order-2 and order-3 sums become sums of zeros, so symmetry
  observations are plain monomials of the mirrored offset;
- the fixed-value boundary observation $u_b(\boldsymbol{x}_q)-u_P$ loses its
  $-M_P^{(\alpha)}$ term;
- `cellValueCoeffsAtPoint` loses the second- and third-derivative moment
  contractions, leaving $u(\boldsymbol{x})=u_P+\sum_\alpha
  \boldsymbol{r}^{\alpha}/\alpha!\;\partial^\alpha u_P$, with the owner
  coefficient still $1-\sum_N C_N$;
- gradients (`cellGradCoeffsAtPoint`, `grad`, `secondGrad`, `thirdGrad`,
  `calcFaceGradCoeffs` owner/neighbour averaging, `fGrad`, the processor-face
  exchange, `faceCentreValues` for the alpha stabilisation) never used moments
  and are unchanged.

The unknowns are the same $N_p$ derivatives; `minNn()`, the stencils, the
weights, the $h$ scaling, the QR solve and every coefficient layout are
reused unchanged.

## 3. Class design

### 3.1 `leastSquaresScheme` (interface)

One addition:

```cpp
        //- True when the stored cell values are cell averages, false when
        //  they are point values at the cell centres
        virtual bool cellAverageUnknowns() const = 0;
```

Results: `movingLeastSquares` false, `leastSquares` false,
`kExactLeastSquares` true. Nothing else in the interface changes.

### 3.2 `cellLeastSquares` (new abstract base)

Created from the present `kExactLeastSquares` implementation with `git mv`,
so the diff shows the edits rather than a 2500-line addition:

- `kExactLeastSquares.H/.C/Templates.C` become
  `cellLeastSquares/cellLeastSquares.H/.C/Templates.C`, class renamed,
  `TypeName` and `addToRunTimeSelectionTable` removed, constructor made
  protected, `cellAverageUnknowns()` left pure virtual.
- It owns everything the two schemes share: stencil, quadrature and weight
  function pointers and their `make*()`; `polynomialOrder_`,
  `cellStencilExtraCells_`, `haloDepthScale_`, `weightFunctionCoeffs_`,
  `includePatchInStencils_`, `calcConditionNumber_`; exponent generation,
  `rowOf`, `derivativeRows`, `minNn`; `QRSolve`; `calcCellCoeffs`,
  `calcFaceGradCoeffs`, `makeFaceGradStencil`; the cell derivative coefficient
  tables and the face-gradient, boundary-data and condition-number tables;
  `cellGradCoeffsAtPoint`, `cellValueCoeffsAtPoint`, `evaluateAtPoint`,
  `valueAtPoint`, `grad`, `secondGrad`, `thirdGrad`, `fGrad`,
  `patchFaceQuadValues`, `exchangeProcessorPatchField`, `clear()`; the two
  derivative-moment contraction statics; the `polynomialOrder()`,
  `quadrature()`, `stencil()`, `faceGradStencil()`, `faceGradCoeffs()`
  accessors.
- Representation switch. The only difference between the two schemes is how
  a central moment of positive order is evaluated. Inside
  `calcCellCoeffs()`, `calcFaceGradCoeffs()`, `cellValueCoeffsAtPoint()` and
  `evaluateAtPoint()` the base calls `cellAverageUnknowns()` once, stores it in
  a local `const bool`, and:
  - the two `centralMoment` lambdas return
    `(i + j + k == 0) ? 1.0 : 0.0` when it is false, before touching any
    table; the moment tables and their remote exchange (`remoteFieldMap`) are
    only built when it is true;
  - `cellValueCoeffsAtPoint` leaves the moment pointers null and skips the two
    contraction terms when it is false;
  - `evaluateAtPoint` requests `quadrature().secondOrderCellMoments()` and
    `thirdOrderCellMoments()` only when it is true.

  The switch is evaluated once per coefficient build, never inside the hot
  loops, so the k-exact path keeps its current cost and stays bit-identical.
  A per-observation virtual hook (`observationMonomial()` or similar) was
  considered and rejected: it would sit in the innermost loop that the pow
  work in #663 just optimised, for two implementations that differ by one
  line.
- No k-exact-specific public API: the moment accessors
  `secondOrderCellMoments()` and `thirdOrderCellMoments()` are not on the
  base; the base reads them from `quadrature()` internally.

### 3.3 `kExactLeastSquares` (thin)

Keeps its directory, `TypeName("kExactLeastSquares")`, its registration, its
README and its `Description`. The class body reduces to: the public
constructor forwarding to `cellLeastSquares`, `cellAverageUnknowns()`
returning true, and the two moment accessors
`secondOrderCellMoments()`/`thirdOrderCellMoments()` kept as thin forwards to
`quadrature()` since they are k-exact concepts (their only external user today
is `Test-cellMoments`, which reads the quadrature directly).

### 3.4 `leastSquares` (new, thin)

`higherOrderHelpers/leastSquares/leastSquares.H/.C`:

```cpp
class leastSquares
:
    public cellLeastSquares
{
public:

    //- Runtime type information
    TypeName("leastSquares");

    //- Construct from components
    leastSquares
    (
        const fvMesh& mesh,
        const boolList& includePatchInStencils,
        const dictionary& dict
    );

    //- Destructor
    virtual ~leastSquares();

    //- The stored cell values are point values at the cell centres
    virtual bool cellAverageUnknowns() const
    {
        return false;
    }
};
```

plus `addToRunTimeSelectionTable(leastSquaresScheme, leastSquares,
dictionary)`. The `Description` says the stored value is the point value at
the cell centre, that one polynomial per cell is built and evaluated at face
quadrature points, and that this is the k-exact formulation with all
positive-order moments set to zero.

### 3.5 `movingLeastSquares`

Unchanged apart from implementing `cellAverageUnknowns()` as false. It is
not moved under `cellLeastSquares`. The small utilities it shares with the
new base (exponent generation, `rowOf`, `derivativeRows`, `QRSolve`) stay
duplicated for now; extracting them into standalone helpers is a possible
later step, not part of this work.

## 4. Representation query outside the reconstruction classes

Seven places currently infer the representation from the type: the
`Test-highOrderGrad` utility (type string), the linear-elastic and
neo-Hookean `manufacturedSolution` function objects
(`isA<kExactLeastSquares>`), and the `cantilever2d`, `sphericalCavity`,
`pressurisedCylinder` and `plateHole` analytical function objects (type
string read from the dictionary). All seven switch to
`scheme.cellAverageUnknowns()`; the `isA` checks and the
`#include "kExactLeastSquares.H"` lines in the two manufactured-solution
function objects go. `Test-highOrderGrad` drops its supported-type check,
since any `leastSquaresScheme` now answers the question.

## 5. Files

New:

- `src/solids4FoamModels/higherOrderHelpers/cellLeastSquares/cellLeastSquares.H`
- `src/solids4FoamModels/higherOrderHelpers/cellLeastSquares/cellLeastSquares.C`
  (via `git mv` from `kExactLeastSquares.C`)
- `src/solids4FoamModels/higherOrderHelpers/cellLeastSquares/cellLeastSquaresTemplates.C`
  (via `git mv`)
- `src/solids4FoamModels/higherOrderHelpers/leastSquares/leastSquares.H`
- `src/solids4FoamModels/higherOrderHelpers/leastSquares/leastSquares.C`
- `src/solids4FoamModels/higherOrderHelpers/leastSquares/README.md`

Modified:

- `kExactLeastSquares.H/.C`: thin class as in 3.3
- `kExactLeastSquares/README.md`: see section 7
- `leastSquaresScheme.H`: `cellAverageUnknowns()`
- `movingLeastSquares.H`: `cellAverageUnknowns()`
- `src/solids4FoamModels/Make/files.openfoam` and `Make/files.foamextend`:
  add `cellLeastSquares.C` and `leastSquares.C` before
  `kExactLeastSquares.C`
- `applications/test/highOrderGrad/Test-highOrderGrad.C`
- the six tutorial function objects of section 4
- regression scripts and tutorial dictionaries of section 6
- `CHANGELOG.md`, `KEXACT_CLAUDE_HANDOFF.md`

## 6. Tests and regressions

- `Test-highOrderGrad` with `type leastSquares;` on the `cooksMembrane`
  (2-D) and `sphericalCavity` (3-D) meshes, serial and two ranks, orders 1 to
  3. Expected: exact reproduction of polynomials up to the reconstruction
  order from point values, as for `movingLeastSquares`.
- After the `git mv` and before adding `leastSquares`: rerun the k-exact
  `plateHole` and `sphericalCavity` regressions and `Test-highOrderGrad` and
  confirm the numbers are unchanged to all printed digits.
- Regression scripts: add `highOrder-leastSquares` to `APPROACHES` in
  `cantilever2d`, `cooksMembrane`, `pressurisedCylinder`, `plateHole`,
  `sphericalCavity`, and `highOrder-leastSquares-parallel` in `plateHole` and
  `sphericalCavity`; extend each `sed` alternation to
  `(movingLeastSquares|kExactLeastSquares|leastSquares)`. Add the approach to
  both `manufacturedSolution` regressions, expecting the "point-valued
  analytical displacement" message. Tolerances are measured first and widened
  only if needed; the expectation is accuracy between `movingLeastSquares`
  and `kExactLeastSquares`.
- CI cost: one more run per case plus two parallel runs, roughly +15% on the
  regression jobs.
- Compile on OpenFOAM v2412, OpenFOAM 9 and foam-extend 4.1.
- Timings on the trigonometric manufactured solution (30k cells, p=3, serial
  and 4 ranks) for the three schemes; `leastSquares` should have the cheapest
  setup since it builds no moments.

## 7. Documentation

- `leastSquares/README.md`: what the scheme is, the point-value statement,
  the zero-moment reduction of section 2, the dictionary entries, and a
  pointer to the k-exact README for the shared machinery.
- `kExactLeastSquares/README.md`: stays where it is. Its function-by-function
  sections (`calcCellCoeffs`, `calcFaceGradCoeffs`, ...) now describe code in
  `cellLeastSquares`; add a short section at the top on the class split and
  say that the derivation is written for cell averages and reduces to the
  point-value scheme as in section 2. Moving the file to
  `cellLeastSquares/README.md` is an alternative (open question 2).
- `CHANGELOG.md`, Unreleased, Added: "Added the `leastSquares` high-order
  reconstruction scheme: point values at the cell centres, one polynomial per
  cell evaluated at the face quadrature points, selected with
  `type leastSquares;`. The per-cell-polynomial machinery it shares with
  `kExactLeastSquares` moved to the new `cellLeastSquares` base, and
  `leastSquaresScheme::cellAverageUnknowns()` reports whether a scheme stores
  cell averages."
- `KEXACT_CLAUDE_HANDOFF.md`: names list, ownership diagram, "important
  source files".

## 8. Order of work

1. Branch from `development`.
2. `cellAverageUnknowns()` on `leastSquaresScheme` and in
   `movingLeastSquares` and `kExactLeastSquares`; switch the seven callers of
   section 4. Compile, run `Test-highOrderGrad` for both existing schemes.
3. `git mv` the k-exact sources to `cellLeastSquares`, rename the class,
   make it abstract, add the representation switch, reduce
   `kExactLeastSquares` to the thin class, update both `Make/files`. Compile
   on the three forks; confirm k-exact results are bit-identical.
4. Add `leastSquares`; compile; `Test-highOrderGrad` serial and parallel.
5. Regression scripts and manufactured-solution regressions; run them and
   set tolerances.
6. READMEs, CHANGELOG, handoff; markdownlint.
7. PR with the timings table.

Steps 2 and 3 are each self-contained and could be separate PRs if the
reviewers prefer smaller diffs.

Estimated size: the `git mv` plus about 60 changed lines in the moved
files, 80 lines for the thin `kExactLeastSquares`, 120 lines for
`leastSquares`, 40 lines in the test and function objects, 60 lines in
regression scripts, plus documentation.

## 9. Open questions

1. Confirm the representation switch of 3.2 (a local `const bool` from
   `cellAverageUnknowns()` inside the base) rather than per-observation
   virtual hooks.
2. Keep the k-exact README in `kExactLeastSquares/` with a class-split note
   (proposed), or `git mv` it to `cellLeastSquares/README.md`.
3. Whether steps 2 and 3 of section 8 should be separate PRs ahead of the new
   scheme.
4. Whether to add `solidProperties.highOrder-leastSquares` files to the two
   `manufacturedSolution` tutorials so their verification sweeps can run it.
