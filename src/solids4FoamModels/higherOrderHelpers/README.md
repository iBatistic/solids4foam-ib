# Higher-order least-squares reconstruction

This directory contains the quadrature, stencil and reconstruction tools used
by the higher-order finite-volume discretisation in solids4foam.

This page records both the current reconstruction framework and the k-exact
scheme under development. The mesh-object manager, abstract scheme interface
and `movingLeastSquares` implementation are in place. `kExactLeastSquares` now
constructs and evaluates cell-centred derivative coefficients, internal- and
processor-face gradients, and weighted boundary gradients on
`fixedDisplacement` and `fixedRotation` boundaries. Traction boundaries are
skipped because their final quadrature
tractions are prescribed directly by the boundary condition. Symmetry and
other fixed-value boundary-condition types are not implemented yet.

## Motivation

The existing `movingLeastSquares` reconstruction treats a cell value as the
value of the reconstructed field at the cell centre. It determines a local
polynomial by solving a weighted least-squares problem at each reconstruction
location.

This interpretation is suitable for the existing pointwise reconstruction,
but a transient finite-volume formulation naturally stores the cell average:

$$
\overline{\phi}_i = \frac{1}{V_i}\int_{V_i}\phi(\boldsymbol{x})\,\mathrm{d}V.
$$

Using the cell average directly as the unknown allows the time derivative to
be formed from stored cell values. Otherwise, the reconstructed field must be
integrated over every cell whenever the time term is evaluated.

When selected, the k-exact reconstruction interprets a `volField` internal
value as a cell average. Both interpretations require similar mesh stencils,
quadrature data and parallel communication, but they use different
reconstruction equations and remain separate schemes.

## Class hierarchy

The implemented hierarchy is:

```text
fvMesh registry
└── leastSquaresReconstruction        concrete mesh-object manager
    ├── displacement scheme           selected at run time
    │   └── movingLeastSquares | kExactLeastSquares
    └── pressure scheme               selected at run time
        └── movingLeastSquares | kExactLeastSquares

leastSquaresScheme                    abstract reconstruction interface
├── movingLeastSquares                pointwise weighted least squares
└── kExactLeastSquares                cell-average k-exact reconstruction

leastSquaresStencil                   shared stencil and parallel exchange
fvMeshQuadrature                      quadrature points and physical weights
```

`leastSquaresReconstruction` and `leastSquaresScheme` have different jobs.
The former is a concrete mesh object that owns reconstruction schemes. The
latter is an abstract run-time selectable interface implemented by the two
numerical methods.

The manager cannot itself be the abstract base class. OpenFOAM's
`MeshObject::New(mesh)` must be able to construct the registered concrete type
when the object is first requested.

## Class responsibilities

### `leastSquaresReconstruction`

`leastSquaresReconstruction` is the single reconstruction object registered
on an `fvMesh`. It:

- own field-specific `leastSquaresScheme` instances;
- create a scheme lazily when it is first requested;
- use the field role, such as `displacement` or `pressure`, as the lookup key;
- retain each scheme for reuse by all models operating on the same mesh;
- forward mesh-motion and topology-change notifications to its schemes.

The current manager contains explicit displacement and pressure pointers. A
general pointer table can be introduced later if arbitrary field roles are
needed.

### `leastSquaresScheme`

`leastSquaresScheme` defines operations required by the solvers and
higher-order operators, independent of the reconstruction method. Its
interface contains only capabilities shared by MLS and k-exact
reconstruction, for example:

```cpp
class leastSquaresScheme
{
public:

    virtual ~leastSquaresScheme() = default;

    virtual label polynomialOrder() const = 0;

    virtual const fvMeshQuadrature& quadrature() const = 0;

    virtual const leastSquaresStencil& stencil() const = 0;

    virtual const CompactListList<label>& faceGradStencil() const = 0;

    virtual const List<CompactListList<vector>>&
    faceGradCoeffs() const = 0;

    virtual void clear() const = 0;
};
```

The interface also contains explicit virtual overloads of `grad()`, `fGrad()`,
`secondGrad()`, `thirdGrad()` and `patchFaceQuadValues()` for the scalar and
vector field types required by the solvers. Templates cannot be virtual, so
each derived class uses templates internally to avoid duplicating its
implementation.

Implementation-specific operations, such as access to the MLS weight
function, remain in the derived class. Solver code normally uses only the
common interface.

`faceGradStencil()` and `faceGradCoeffs()` form one addressing-coefficient
pair. `movingLeastSquares` returns its existing face stencil.
`kExactLeastSquares` combines the owner- and neighbour-cell reconstruction
stencils on internal faces and uses a one-sided owner stencil on prescribed
displacement boundaries. Both `fGrad()` and `hofvm` use this pair, ensuring
that the explicit
residual and implicit Jacobian apply the same cell-unknown face-gradient
operator.

### `movingLeastSquares`

`movingLeastSquares` remains responsible for the existing pointwise method:

- construction of the weighted least-squares systems;
- QR solution and optional condition-number calculation;
- pointwise interpolation and derivative coefficients;
- MLS weight-function selection;
- MLS-specific cached coefficient data.

It derives from `leastSquaresScheme` without changing its numerical method.

### `kExactLeastSquares`

`kExactLeastSquares` is responsible for the cell-average method:

- construction of a polynomial whose cell averages satisfy the k-exact
  constraints;
- access to volume-normalised central cell moments owned and calculated by
  `fvMeshQuadrature`;
- cell derivative coefficients and gradients at face quadrature points;
- weighted `fixedDisplacement` and `fixedRotation` boundary gradients;
- k-exact-specific conditioning and consistency checks.

The central cell average is imposed by the k-exact formulation. Neighbouring
cell averages provide the remaining least-squares equations.

### `leastSquaresStencil`

`leastSquaresStencil` contains reconstruction-independent mesh operations:

- construction of cell and face stencils;
- global cell indexing;
- discovery and storage of remote stencil cells;
- parallel exchange of remote field values and cell centres;
- clearing stencil data after a mesh change.

The current cell stencil excludes the central cell. The MLS formulation adds
the central-cell coefficient separately. The k-exact formulation instead
multiplies every stored neighbour coefficient by the difference between the
neighbour and central cell averages.

### `fvMeshQuadrature`

Each concrete reconstruction scheme owns its own `fvMeshQuadrature` instance:

```text
movingLeastSquares
├── leastSquaresStencil
└── fvMeshQuadrature

kExactLeastSquares
├── leastSquaresStencil
└── fvMeshQuadrature
```

Separate ownership allows displacement and pressure schemes to use different
polynomial and integration orders. It also avoids coupling the two numerical
methods through a shared mutable cache.

`fvMeshQuadrature` calculates and caches quadrature points, their physical area
or volume weights, and raw volume-normalised central cell moments. First-order
moments are created for integration order one or higher, second-order moments
for order two or higher, and third-order moments only for order three. In two
dimensions, the empty direction is omitted and its moment components remain
exactly zero.

In parallel, the two copies of a processor face can generate the same
quadrature points in a different order because their face-vertex order is
opposite. The lower-rank processor is therefore the master for each processor
patch. On the first call to `faceQuadPoints()` or `faceQuadWeights()`, it sends
the number of points per face followed by its point and weight arrays. The
higher-rank processor checks that every corresponding face has the same number
of points and directly replaces its local arrays with the master's arrays.
Thus, a quadrature-point index identifies the same physical point and weight
on both sides of every processor face. This synchronisation is performed once
and reused with the cached quadrature data.

Scaled, shifted or stencil-exchanged moment data remain the responsibility of
`kExactLeastSquares`, because those representations depend on the selected
polynomial basis and reconstruction procedure.

`fvMeshQuadrature` does not need to become another mesh object. It is owned by
a concrete scheme, and its cached geometry is cleared through that scheme when
the mesh changes.

## K-exact cell coefficient construction

This section describes the equations implemented by
`kExactLeastSquares::calcCellCoeffs()`. It is intended to allow a direct
comparison between the code and the mathematical derivation.

### Read this first: the derivation in explicit notation

This subsection follows the notation in **eq_derivation.tex**: ordinary
polynomial coefficients are \(c_x,c_{xx},\ldots\). The compact notation used
later is useful in generic loops, but it hides the individual matrix entries.
We therefore start with the complete two-dimensional quadratic equation.

#### Step 1: eliminate the constant polynomial coefficient

For owner cell \(P\), define

$$
r_x=x-x_P,\qquad r_y=y-y_P.
$$

An ordinary quadratic reconstruction is

$$
u_P^R(\boldsymbol{x})
=c_0+c_xr_x+c_yr_y+c_{xx}r_x^2+c_{xy}r_xr_y+c_{yy}r_y^2.
$$

The stored value \(\bar{u}_P\) is the cell average, not
\(u(\boldsymbol{x}_P)\). Define

$$
M_{xx}^{P}
=\frac{1}{|\Omega_P|}\int_{\Omega_P}(x-x_P)^2\,\mathrm{d}V,
$$

$$
M_{xy}^{P}
=\frac{1}{|\Omega_P|}\int_{\Omega_P}
(x-x_P)(y-y_P)\,\mathrm{d}V,
$$

$$
M_{yy}^{P}
=\frac{1}{|\Omega_P|}\int_{\Omega_P}(y-y_P)^2\,\mathrm{d}V.
$$

The first moments vanish because \(\boldsymbol{x}_P\) is the volume centroid.
Averaging the ordinary polynomial over \(P\) gives

$$
\bar{u}_P
=c_0+c_{xx}M_{xx}^{P}+c_{xy}M_{xy}^{P}+c_{yy}M_{yy}^{P}.
$$

Therefore,

$$
c_0
=\bar{u}_P-c_{xx}M_{xx}^{P}-c_{xy}M_{xy}^{P}-c_{yy}M_{yy}^{P}.
$$

Substituting this expression back into the polynomial gives

$$
\begin{aligned}
u_P^R(\boldsymbol{x})={}&
\bar{u}_P+c_xr_x+c_yr_y\\
&+c_{xx}\left(r_x^2-M_{xx}^{P}\right)\\
&+c_{xy}\left(r_xr_y-M_{xy}^{P}\right)\\
&+c_{yy}\left(r_y^2-M_{yy}^{P}\right).
\end{aligned}
$$

This is the essential k-exact step. It removes \(c_0\) from the unknowns and
guarantees that the average of the reconstruction over \(P\) is always
\(\bar{u}_P\), whatever values the remaining coefficients take.

#### Step 2: average the owner reconstruction over neighbour N

For one stencil neighbour \(N\), define

$$
d_x=x_N-x_P,\qquad d_y=y_N-y_P.
$$

Inside \(N\),

$$
x-x_P=(x-x_N)+d_x,\qquad y-y_P=(y-y_N)+d_y.
$$

The average of a linear term is simply

$$
\frac{1}{|\Omega_N|}\int_{\Omega_N}(x-x_P)\,\mathrm{d}V=d_x,
$$

because the average of \(x-x_N\) is zero. For the \(xx\) term,

$$
(x-x_P)^2=(x-x_N)^2+2d_x(x-x_N)+d_x^2.
$$

The middle term has zero average, so

$$
\frac{1}{|\Omega_N|}\int_{\Omega_N}(x-x_P)^2\,\mathrm{d}V
=d_x^2+M_{xx}^{N}.
$$

Similarly,

$$
\frac{1}{|\Omega_N|}\int_{\Omega_N}
(x-x_P)(y-y_P)\,\mathrm{d}V=d_xd_y+M_{xy}^{N},
$$

and

$$
\frac{1}{|\Omega_N|}\int_{\Omega_N}(y-y_P)^2\,\mathrm{d}V
=d_y^2+M_{yy}^{N}.
$$

The equation contributed by neighbour \(N\) is therefore

$$
\begin{aligned}
\bar{u}_N-\bar{u}_P={}&c_xd_x+c_yd_y\\
&+c_{xx}\left(d_x^2+M_{xx}^{N}-M_{xx}^{P}\right)\\
&+c_{xy}\left(d_xd_y+M_{xy}^{N}-M_{xy}^{P}\right)\\
&+c_{yy}\left(d_y^2+M_{yy}^{N}-M_{yy}^{P}\right).
\end{aligned}
$$

The practical rule for checking a quadratic entry is

$$
\boxed{
\text{quadratic entry}
=\text{centre-offset term}
+\text{moment of }N
-\text{moment of }P
}.
$$

#### Step 3: form the unscaled matrix

Using

$$
\boldsymbol{c}_P=
\begin{bmatrix}c_x&c_y&c_{xx}&c_{xy}&c_{yy}\end{bmatrix}^{T},
$$

one neighbour contributes the row

$$
\boldsymbol{G}_N=
\begin{bmatrix}
d_x,\;
d_y,\;
d_x^2+M_{xx}^{N}-M_{xx}^{P},\;
d_xd_y+M_{xy}^{N}-M_{xy}^{P},\;
d_y^2+M_{yy}^{N}-M_{yy}^{P}
\end{bmatrix}.
$$

Thus,

$$
\boldsymbol{G}_N\boldsymbol{c}_P=\bar{u}_N-\bar{u}_P.
$$

For \(m\) neighbours,

$$
\underbrace{
\begin{bmatrix}
\boldsymbol{G}_{N_1}\\
\boldsymbol{G}_{N_2}\\
\vdots\\
\boldsymbol{G}_{N_m}
\end{bmatrix}}_{\boldsymbol{G}_P\;(m\times5)}
\boldsymbol{c}_P
=
\underbrace{
\begin{bmatrix}
\bar{u}_{N_1}-\bar{u}_P\\
\bar{u}_{N_2}-\bar{u}_P\\
\vdots\\
\bar{u}_{N_m}-\bar{u}_P
\end{bmatrix}}_{\boldsymbol{b}_P\;(m\times1)}.
$$

This is the matrix called \(\boldsymbol{A}_P\) in
**eq_derivation.tex**. We call it \(\boldsymbol{G}_P\) here because the C++
code uses the name **A** for the final reconstruction matrix.

#### Step 4: add the cubic columns

For \(p=3\), append the mean-free terms

$$
\begin{aligned}
&c_{xxx}(r_x^3-M_{xxx}^{P})
+c_{xxy}(r_x^2r_y-M_{xxy}^{P})\\
&\qquad
+c_{xyy}(r_xr_y^2-M_{xyy}^{P})
+c_{yyy}(r_y^3-M_{yyy}^{P})
\end{aligned}
$$

to the polynomial. The four additional entries in
\(\boldsymbol{G}_N\) are

$$
G_{N,xxx}
=d_x^3+3d_xM_{xx}^{N}+M_{xxx}^{N}-M_{xxx}^{P},
$$

$$
G_{N,xxy}
=d_x^2d_y+d_yM_{xx}^{N}+2d_xM_{xy}^{N}
+M_{xxy}^{N}-M_{xxy}^{P},
$$

$$
G_{N,xyy}
=d_xd_y^2+d_xM_{yy}^{N}+2d_yM_{xy}^{N}
+M_{xyy}^{N}-M_{xyy}^{P},
$$

$$
G_{N,yyy}
=d_y^3+3d_yM_{yy}^{N}+M_{yyy}^{N}-M_{yyy}^{P}.
$$

For example,

$$
\begin{aligned}
(r_{N,x}+d_x)^2(r_{N,y}+d_y)={}&
r_{N,x}^2r_{N,y}+d_y r_{N,x}^2
+2d_xr_{N,x}r_{N,y}\\
&+2d_xd_y r_{N,x}+d_x^2r_{N,y}+d_x^2d_y.
\end{aligned}
$$

After averaging over \(N\), the two first-moment terms vanish and the
remaining terms give \(G_{N,xxy}\). The nested loops in
**averageMonomial()** perform these same expansions for all exponents.

#### Why averageMonomial() contains three nested loops

The three loops do not represent three cells or three polynomial orders. They
represent the three coordinate directions \(x\), \(y\), and \(z\).

Suppose the requested owner-centred monomial is

$$
(x-x_P)^i(y-y_P)^j(z-z_P)^k.
$$

Within cell \(N\), each coordinate is split into a coordinate relative to the
centre of \(N\) and a fixed centre offset:

$$
x-x_P=r_{N,x}+d_x,
\qquad
y-y_P=r_{N,y}+d_y,
\qquad
z-z_P=r_{N,z}+d_z.
$$

Consequently, the required monomial is a product of three binomials:

$$
(r_{N,x}+d_x)^i
(r_{N,y}+d_y)^j
(r_{N,z}+d_z)^k.
$$

Each binomial has its own summation index:

- \(a=0,\ldots,i\) expands the \(x\) factor;
- \(b=0,\ldots,j\) expands the \(y\) factor;
- \(c=0,\ldots,k\) expands the \(z\) factor.

Expanding all three factors gives

$$
\begin{aligned}
&(r_{N,x}+d_x)^i
(r_{N,y}+d_y)^j
(r_{N,z}+d_z)^k\\
&\quad =
\sum_{a=0}^{i}
\sum_{b=0}^{j}
\sum_{c=0}^{k}
\binom{i}{a}\binom{j}{b}\binom{k}{c}
d_x^{i-a}d_y^{j-b}d_z^{k-c}
r_{N,x}^{a}r_{N,y}^{b}r_{N,z}^{c}.
\end{aligned}
$$

When this expression is averaged over cell \(N\), the last product becomes
the central moment

$$
M_{abc}^{N}
=
\frac{1}{|\Omega_N|}
\int_{\Omega_N}
r_{N,x}^{a}r_{N,y}^{b}r_{N,z}^{c}\,\mathrm{d}V.
$$

Therefore, **averageMonomial()** calculates

$$
\sum_{a=0}^{i}
\sum_{b=0}^{j}
\sum_{c=0}^{k}
\binom{i}{a}\binom{j}{b}\binom{k}{c}
d_x^{i-a}d_y^{j-b}d_z^{k-c}M_{abc}^{N}.
$$

This equation maps directly to the three loops:

| Loop variable | Coordinate | Selects |
|---|---|---|
| \(a\) | \(x\) | power of \(r_{N,x}\) placed in the cell moment |
| \(b\) | \(y\) | power of \(r_{N,y}\) placed in the cell moment |
| \(c\) | \(z\) | power of \(r_{N,z}\) placed in the cell moment |

The remaining powers, \(i-a\), \(j-b\), and \(k-c\), are powers of the
centre offsets \(d_x\), \(d_y\), and \(d_z\).

For example, the **xxy** monomial has

$$
i=2,\qquad j=1,\qquad k=0.
$$

The loop ranges are then

$$
a=0,1,2,\qquad b=0,1,\qquad c=0.
$$

There are \(3\times2\times1=6\) combinations. They generate the six terms in

$$
(r_{N,x}+d_x)^2(r_{N,y}+d_y).
$$

Terms containing \(M_{100}^{N}\) or \(M_{010}^{N}\) later vanish because first
central moments are zero, but it is simpler and less error-prone for the
generic loop to generate them and let **centralMoment()** return zero.

For a two-dimensional mesh, **calcCellCoeffs()** sets \(k=0\). The third loop
then executes only once, with \(c=0\), and contributes the factor
\(\binom{0}{0}d_z^0=1\). Thus, the same implementation handles both 2-D and
3-D without separate expansion code.

When **averageMonomial()** is called for the owner cell, the offset
\(\boldsymbol{d}\) is zero. All terms containing a positive power of an
offset vanish, leaving only \(a=i\), \(b=j\), and \(c=k\). The result is then
exactly the owner central moment \(M_{ijk}^{P}\).

#### Step 5: connect the ordinary coefficients to code derivatives

The derivation uses ordinary polynomial coefficients. At the owner centre,

$$
c_x=u_x,\qquad c_y=u_y,
$$

$$
c_{xx}=\frac{u_{xx}}{2},\qquad
c_{xy}=u_{xy},\qquad
c_{yy}=\frac{u_{yy}}{2},
$$

and

$$
c_{xxx}=\frac{u_{xxx}}{6},\qquad
c_{xxy}=\frac{u_{xxy}}{2},\qquad
c_{xyy}=\frac{u_{xyy}}{2},\qquad
c_{yyy}=\frac{u_{yyy}}{6}.
$$

The implementation solves for derivatives, not for the ordinary \(c\)
coefficients. It also scales each derivative using

$$
h_P=2\max_{N\in\mathcal{S}_P}
|\boldsymbol{x}_N-\boldsymbol{x}_P|.
$$

For example, the unknown corresponding to \(xx\) is \(h_P^2u_{xx}\).
Its matrix entry is therefore

$$
Q_{xx,N}
=\frac{G_{N,xx}}{2h_P^2}
=\frac{d_x^2+M_{xx}^{N}-M_{xx}^{P}}{2h_P^2}.
$$

The product is unchanged:

$$
\frac{G_{N,xx}}{2h_P^2}(h_P^2u_{xx})
=G_{N,xx}\frac{u_{xx}}{2}
=G_{N,xx}c_{xx}.
$$

For the mixed quadratic term there is no factor two:

$$
Q_{xy,N}
=\frac{d_xd_y+M_{xy}^{N}-M_{xy}^{P}}{h_P^2},
$$

because \(c_{xy}=u_{xy}\). For a cubic mixed term,

$$
Q_{xxy,N}=\frac{G_{N,xxy}}{2h_P^3},
$$

because \(c_{xxy}=u_{xxy}/2\). The factorials do not change the
reconstruction; they make the solved coefficients derivatives directly.

#### Step 6: understand the matrix orientation in calcCellCoeffs()

The natural scaled equation has one row per neighbour:

$$
\boldsymbol{G}^{\mathrm{scaled}}_P\boldsymbol{a}_P
=\boldsymbol{b}_P.
$$

The C++ variable **Q** stores its transpose:

$$
\mathtt{Q}
=\left(\boldsymbol{G}^{\mathrm{scaled}}_P\right)^T.
$$

Consequently:

- **Q.rows()** is the number of non-constant polynomial terms;
- **Q.cols()** is the number of stencil neighbours;
- **Q(p, cI)** is non-constant term \(p\) for neighbour **cI**;
- no row offset is needed because the exponent list excludes the constant term.

The number of reconstructed unknowns is:

| Order | Two dimensions | Three dimensions |
|------:|---------------:|-----------------:|
| 1 | 2 | 3 |
| 2 | 5 | 9 |
| 3 | 9 | 19 |

For k-exact reconstruction, **minNn()** is the number of non-constant
polynomial terms shown in the table. This is also the minimum number of
neighbour equations needed by the least-squares system. The
**leastSquaresStencil** constructor uses a different convention: its requested
cell-stencil size includes the central cell. Therefore **makeStencils()** adds
one only when passing the size to **leastSquaresStencil**. The returned
**cellsStencil()** contains **minNn() + cellStencilExtraCells** neighbours.

#### Step 7: understand the weighting convention

Let

$$
\boldsymbol{W}
=\operatorname{diag}(w_{N_1},\ldots,w_{N_m}),
$$

where **weightFunc().weight()** supplies \(w_N\). The code minimizes

$$
\left\|
\boldsymbol{W}^{1/2}
\left(
\boldsymbol{G}^{\mathrm{scaled}}_P\boldsymbol{a}
-\boldsymbol{b}_P
\right)
\right\|_2^2,
$$

or, equivalently,

$$
\sum_{N\in\mathcal{S}_P}
w_N
\left(
\boldsymbol{G}^{\mathrm{scaled}}_N\boldsymbol{a}
-(\bar{u}_N-\bar{u}_P)
\right)^2.
$$

This explains the **cwiseSqrt()** call in **QRSolve()**. In
**eq_derivation.tex**, the matrix multiplying the residual is denoted by
\(\boldsymbol{W}_P\). The notations agree if that matrix is interpreted as
the square root of the weight matrix above. Otherwise,
\(\|\boldsymbol{W}r\|^2\) squares each numerical weight a second time.

**QRSolve()** returns the geometry-only matrix **A**:

$$
\boldsymbol{a}_P=\mathtt{A}\boldsymbol{b}_P.
$$

Thus, **A** in the code corresponds to \(\boldsymbol{C}_P\) in
**eq_derivation.tex**; it is not the neighbour-row matrix called
\(\boldsymbol{A}_P\) in that document.

#### Step 8: follow A into the derivative evaluations

Rows of **A** return scaled derivatives. Before storing them,
**calcCellCoeffs()** divides:

- first-derivative rows by \(h_P\);
- second-derivative rows by \(h_P^2\);
- third-derivative rows by \(h_P^3\).

For example,

$$
u_x=\sum_N\frac{A_{xN}}{h_P}(\bar{u}_N-\bar{u}_P),
$$

$$
u_{xx}=\sum_N\frac{A_{xx,N}}{h_P^2}(\bar{u}_N-\bar{u}_P),
$$

and

$$
u_{xxy}=\sum_N\frac{A_{xxy,N}}{h_P^3}
(\bar{u}_N-\bar{u}_P).
$$

Each coefficient list has exactly one entry per stencil neighbour. There is
no stored owner coefficient. The owner contribution appears when **grad()**,
**secondGrad()**, or **thirdGrad()** forms
\(\bar{u}_N-\bar{u}_P\).

These functions currently return derivatives at the owner cell centre.
Face-quadrature evaluation differentiates the complete polynomial at the
quadrature point; it does not solve another least-squares problem there.

#### Gradient at an internal-face quadrature point

Consider an internal face \(f\). The mesh addressing identifies its two sides:

$$
P=\mathtt{mesh.owner()[faceI]},
\qquad
N=\mathtt{mesh.neighbour()[faceI]}.
$$

Let \(\boldsymbol{x}_q\) be one quadrature point on the face, and define

$$
\boldsymbol{r}_{Pq}=\boldsymbol{x}_q-\boldsymbol{x}_P,
\qquad
\boldsymbol{r}_{Nq}=\boldsymbol{x}_q-\boldsymbol{x}_N.
$$

The gradient of the reconstruction associated with \(P\) is evaluated at
\(\boldsymbol{x}_q\). For polynomial order one,

$$
\nabla u_P^R(\boldsymbol{x}_q)=\boldsymbol{g}_P,
$$

where \(\boldsymbol{g}_P=\nabla u_P^R(\boldsymbol{x}_P)\). For order two,

$$
\nabla u_P^R(\boldsymbol{x}_q)
=
\boldsymbol{g}_P
+\boldsymbol{H}_P\boldsymbol{r}_{Pq},
$$

where \(\boldsymbol{H}_P\) is the Hessian at the cell centre. For order three,

$$
\nabla u_P^R(\boldsymbol{x}_q)
=
\boldsymbol{g}_P
+\boldsymbol{H}_P\boldsymbol{r}_{Pq}
+\frac{1}{2}
\boldsymbol{T}_P:
\left(
\boldsymbol{r}_{Pq}\otimes\boldsymbol{r}_{Pq}
\right),
$$

where \(\boldsymbol{T}_P\) contains the third derivatives. Cell moments do not
appear explicitly in these gradient equations because they are constants with
respect to \(\boldsymbol{x}\), and their derivatives are zero.

In two dimensions, the cubic formula is easier to check component-by-component:

$$
\begin{aligned}
\left.\frac{\partial u_P^R}{\partial x}\right|_{\boldsymbol{x}_q}
={}&u_x
+u_{xx}r_x
+u_{xy}r_y\\
&+\frac{1}{2}u_{xxx}r_x^2
+u_{xxy}r_xr_y
+\frac{1}{2}u_{xyy}r_y^2,
\end{aligned}
$$

and

$$
\begin{aligned}
\left.\frac{\partial u_P^R}{\partial y}\right|_{\boldsymbol{x}_q}
={}&u_y
+u_{xy}r_x
+u_{yy}r_y\\
&+\frac{1}{2}u_{xxy}r_x^2
+u_{xyy}r_xr_y
+\frac{1}{2}u_{yyy}r_y^2,
\end{aligned}
$$

where \(r_x=x_q-x_P\) and \(r_y=y_q-y_P\). The same equations are evaluated
for the reconstruction associated with \(N\), using
\(\boldsymbol{r}_{Nq}\) and the derivatives reconstructed in cell \(N\).

The initial internal-face operator will use the central value

$$
\nabla u_f(\boldsymbol{x}_q)
=
\frac{1}{2}
\nabla u_P^R(\boldsymbol{x}_q)
+\frac{1}{2}
\nabla u_N^R(\boldsymbol{x}_q).
$$

For displacement, these scalar equations are applied to every displacement
component, producing the displacement-gradient tensor.

#### cellGradCoeffsAtPoint

The private helper has the following interface:

```cpp
void cellGradCoeffsAtPoint
(
    const label cellI,
    const point& x,
    UList<vector>& coeffs
) const;
```

It converts the derivative coefficient operators already stored for one cell
into gradient coefficients at an arbitrary point \(\boldsymbol{x}\).

The field values used with these coefficients must be interpreted as cell
averages,

$$
\overline{\boldsymbol{D}}_P
=
\frac{1}{V_P}
\int_{V_P}\boldsymbol{D}\,\mathrm{d}V,
$$

and not as point values sampled at the cell centres.

The output list has the same size and ordering as
**stencil().cellsStencil()[cellI]**. It does not contain an additional owner
entry. For every stencil position \(s\), it calculates
$$
\boldsymbol{C}_{s}^{P}(\boldsymbol{x})
=
\boldsymbol{G}_{s}^{P}
+\boldsymbol{H}_{s}^{P}\boldsymbol{r}
+\frac{1}{2}
\boldsymbol{T}_{s}^{P}:
(\boldsymbol{r}\otimes\boldsymbol{r}),
\qquad
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P.
$$

This collapsed coefficient is obtained by evaluating the gradient of the
owner-cell reconstruction at \(\boldsymbol{x}\):
$$
\nabla\boldsymbol{D}_P^R(\boldsymbol{x})
=
\nabla\boldsymbol{D}_P
+\nabla^2\boldsymbol{D}_P\cdot\boldsymbol{r}
+\frac{1}{2}
\nabla^3\boldsymbol{D}_P:
(\boldsymbol{r}\otimes\boldsymbol{r}).
$$

Consequently:

- for \(p=1\), the reconstructed gradient is constant inside the cell;
- for \(p=2\), the second derivative makes the gradient vary linearly with
  \(\boldsymbol{r}\);
- for \(p=3\), the third derivative adds a quadratic variation of the
  gradient.

The \(p=2\) contribution can therefore be viewed as linear extrapolation of
the cell-centred gradient. The \(p=3\) contribution is a quadratic correction,
so the complete third-order evaluation is not merely linearised
extrapolation.

Here:

- \(\boldsymbol{G}_{s}^{P}\) comes from **cellGradCoeffs()**;
- \(\boldsymbol{H}_{s}^{P}\) comes from **cellSecondGradCoeffs()**;
- \(\boldsymbol{T}_{s}^{P}\) comes from **cellThirdGradCoeffs()**.

The Hessian term is added only for \(p\ge2\), and the third-derivative term
only for \(p\ge3\). In the C++ implementation, the two contractions are

```cpp
secondGradCoeff & r
0.5*((thirdGradCoeff & r) & r)
```

The resulting coefficients retain the neighbour-minus-owner form:

$$
\nabla u_P^R(\boldsymbol{x})
=
\sum_s
\boldsymbol{C}_{s}^{P}(\boldsymbol{x})
\left(
\bar{u}_{J_s}-\bar{u}_P
\right).
$$

For displacement, the same operation is applied component-wise and produces
the displacement-gradient tensor:
$$
\nabla\boldsymbol{D}_P^R(\boldsymbol{x})
=
\sum_s
\boldsymbol{C}_{s}^{P}(\boldsymbol{x})
\otimes
\left(
\overline{\boldsymbol{D}}_{J_s}
-\overline{\boldsymbol{D}}_P
\right).
$$

Thus, multiplying the returned coefficients by cell-average displacement
differences gives the gradient at \(\boldsymbol{x}\). Multiplying them by only
the stencil displacements, without subtracting the owner displacement, would
be incorrect because **cellGradCoeffsAtPoint()** does not return an owner
coefficient.

The helper is independent of the reconstructed field. It combines
geometry-dependent coefficient operators only and can therefore be called
once while constructing the cached face coefficients.

When this difference form is inserted into an ordinary global-cell coefficient
map, every neighbour receives
\(+\boldsymbol{C}_{s}^{P}(\boldsymbol{x})\), while the owner receives

$$
\boldsymbol{C}_{P}^{P}(\boldsymbol{x})
=
-\sum_s\boldsymbol{C}_{s}^{P}(\boldsymbol{x}).
$$

The owner entry is added during face-coefficient construction, not by
**cellGradCoeffsAtPoint()**. After this conversion, the reconstruction can be
applied in the ordinary form

$$
\nabla\boldsymbol{D}_P^R(\boldsymbol{x})
=
\sum_{J\in\{P\}\cup\mathcal{S}_P}
\boldsymbol{C}_{J}^{P}(\boldsymbol{x})
\otimes\overline{\boldsymbol{D}}_J.
$$

Although the helper accepts any point, it is intended for points inside or
close to the owner cell, particularly face quadrature points. Polynomial
extrapolation far outside the reconstruction stencil can be poorly behaved.

#### How owner and neighbour contributions use one face stencil

The owner-side reconstruction can be written as a linear combination of cell
averages:

$$
\nabla u_P^R(\boldsymbol{x}_q)
=
\sum_{J\in\{P\}\cup\mathcal{S}_P}
\boldsymbol{C}_{J}^{P,q}\bar{u}_J.
$$

Similarly,

$$
\nabla u_N^R(\boldsymbol{x}_q)
=
\sum_{J\in\{N\}\cup\mathcal{S}_N}
\boldsymbol{C}_{J}^{N,q}\bar{u}_J.
$$

During coefficient construction, the two sides remain known from
**mesh.owner()** and **mesh.neighbour()**. Temporary owner-side and
neighbour-side coefficient maps can therefore be formed independently.

The stored face stencil is then the union

$$
\mathcal{F}_f
=
\{P\}\cup\mathcal{S}_P
\cup\{N\}\cup\mathcal{S}_N.
$$

The current implementation constructs this union for internal faces in
**makeFaceGradStencil()**. It inserts the global IDs of \(P\) and \(N\), adds
both cell stencils, removes duplicate global IDs, and sorts the result.
Boundary-face rows are intentionally left empty until their boundary
reconstruction is implemented.

#### calcFaceGradCoeffs

For every face quadrature point, **calcFaceGradCoeffs()** converts the two
neighbour-minus-owner reconstructions into one operator that multiplies cell
averages directly.

The coefficients returned by **cellGradCoeffsAtPoint()** for the owner side
satisfy

$$
\nabla u_P^R(\boldsymbol{x}_q)
=
\sum_{J\in\mathcal{S}_P}
\boldsymbol{C}_J^{P,q}
\left(\bar{u}_J-\bar{u}_P\right).
$$

Expanding the differences gives

$$
\begin{aligned}
\nabla u_P^R(\boldsymbol{x}_q)
={}&
\sum_{J\in\mathcal{S}_P}
\boldsymbol{C}_J^{P,q}\bar{u}_J
\\
&-
\left(
\sum_{J\in\mathcal{S}_P}
\boldsymbol{C}_J^{P,q}
\right)\bar{u}_P.
\end{aligned}
$$

Therefore, the explicit coefficient of the owner value is

$$
\boldsymbol{C}_P^{P,q}
=
-\sum_{J\in\mathcal{S}_P}\boldsymbol{C}_J^{P,q}.
$$

The code first accumulates the positive stencil-cell coefficients and their
sum:

```cpp
coeffs[faceStencilI] += 0.5*coeff;
ownSum += coeff;
```

It then inserts the negative owner coefficient:

```cpp
coeffs[ownIndex] -= 0.5*ownSum;
```

Exactly the same expansion is performed for neighbour cell \(N\):

$$
\boldsymbol{C}_N^{N,q}
=
-\sum_{J\in\mathcal{S}_N}\boldsymbol{C}_J^{N,q},
$$

which gives

```cpp
coeffs[neiIndex] -= 0.5*neiSum;
```

The factor \(1/2\) comes from the selected central face gradient

$$
\nabla u_f(\boldsymbol{x}_q)
=
\frac{1}{2}\nabla u_P^R(\boldsymbol{x}_q)
+\frac{1}{2}\nabla u_N^R(\boldsymbol{x}_q).
$$

The code uses `+=` and `-=` instead of assignment because one merged-stencil
entry can receive contributions from both reconstructions. For example, when
\(N\in\mathcal{S}_P\), the final coefficient multiplying \(\bar{u}_N\) is

$$
\boldsymbol{C}_N^{f,q}
=
\frac{1}{2}\boldsymbol{C}_N^{P,q}
-\frac{1}{2}
\sum_{J\in\mathcal{S}_N}\boldsymbol{C}_J^{N,q}.
$$

Similarly, when \(P\in\mathcal{S}_N\),

$$
\boldsymbol{C}_P^{f,q}
=
-\frac{1}{2}
\sum_{J\in\mathcal{S}_P}\boldsymbol{C}_J^{P,q}
+\frac{1}{2}\boldsymbol{C}_P^{N,q}.
$$

If a cell occurs on only one side, the missing contribution from the other
side is zero. If it occurs in both stencils, both contributions are accumulated
at the same global-cell position.

The negative owner and neighbour sums are also necessary for constant-field
preservation. Each one-sided reconstruction satisfies

$$
-\sum_{J\in\mathcal{S}_P}\boldsymbol{C}_J^{P,q}
+\sum_{J\in\mathcal{S}_P}\boldsymbol{C}_J^{P,q}
=\boldsymbol{0},
$$

and similarly for cell \(N\). Consequently, the final merged coefficients
satisfy

$$
\sum_{J\in\mathcal{F}_f}\boldsymbol{C}_J^{f,q}
=\boldsymbol{0},
$$

so a constant displacement field produces zero face gradient. Without
`coeffs[ownIndex] -= 0.5*ownSum` and
`coeffs[neiIndex] -= 0.5*neiSum`, this property would be lost.

The scalar and vector **fGrad()** overloads apply this merged operator directly
to local or remotely collected cell-average values. Internal faces use the
merged operator described here, processor faces exchange and average two
one-sided gradients, and fixed-displacement faces use the weighted one-sided
reconstruction described below. All result rows are first set to zero, so
unsupported physical boundaries remain zero.

For each global cell \(J\) in this union, the final central-gradient
coefficient is

$$
\boldsymbol{C}_{J}^{f,q}
=
\frac{1}{2}\boldsymbol{C}_{J}^{P,q}
+\frac{1}{2}\boldsymbol{C}_{J}^{N,q}.
$$

If \(J\) does not occur in one side's reconstruction, that side contributes
zero. If it occurs in both reconstructions, the two contributions are added
into the same coefficient. The final operator is

$$
\nabla u_f(\boldsymbol{x}_q)
=
\sum_{J\in\mathcal{F}_f}
\boldsymbol{C}_{J}^{f,q}\bar{u}_J.
$$

Consequently, the stored arrays do not require an owner-side or neighbour-side
flag:

- **faceGradStencil()** at indices \((faceI,sI)\) stores the global cell ID
  \(J\);
- **faceGradCoeffs()** at indices \((faceI,qI,sI)\) stores
  \(\boldsymbol{C}_{J}^{f,q}\);
- the common index **sI** associates the global cell with its coefficient.

The side information is required only while constructing the coefficients.
After the two linear reconstructions have been combined, both **fGrad()** and
**hofvm** need only the final global-cell addressing and the final coefficient.
Keeping one merged operator also guarantees that the explicit residual and
implicit Jacobian use exactly the same face gradient.

For a constant field, the final coefficients must satisfy

$$
\sum_{J\in\mathcal{F}_f}\boldsymbol{C}_{J}^{f,q}
=\boldsymbol{0}.
$$

**Test-highOrderGrad** checks exact internal-face gradient reproduction for
linear, quadratic and cubic fields at every internal-face quadrature point.

#### Processor-face gradients

A processor face has one local owner cell on each processor. At a quadrature
point \(\boldsymbol{x}_q\), each processor first evaluates its owner-side
reconstruction using **cellGradCoeffsAtPoint()**:

$$
\boldsymbol{g}_{P,q}
=
\sum_{J\in\mathcal{S}_P}
\boldsymbol{C}^{P,q}_J
(\overline{u}_J-\overline{u}_P).
$$

The two processors exchange these owner-side gradients. The final gradient on
both copies of the processor face is

$$
\boldsymbol{g}_{f,q}
=
\frac{1}{2}
\left(
\boldsymbol{g}_{P,q}+\boldsymbol{g}_{N,q}
\right).
$$

Before either reconstruction is evaluated, `fvMeshQuadrature` gives both
copies of the processor face the same quadrature-point coordinates, weights
and ordering. The lower-rank processor defines this canonical data and the
higher-rank processor stores a copy. Therefore, `fGrad()` exchanges only the
flattened owner-side gradient arrays and pairs entries with the same
quadrature-point index. It does not exchange point counts, coordinates or
weights. Unlike an internal face, a processor face does not store one merged
owner-neighbour coefficient array because the two cell reconstructions are
owned by different processes.

**Test-highOrderGrad** checks both scalar and vector gradients directly at every
processor-face quadrature point.

#### Solid-traction boundary treatment

`solidTraction` prescribes traction, not displacement. The displacement
gradient calculated on that patch is therefore not used to determine the final
traction. The solid model calls `enforceTractionBoundaries()` and replaces the
temporary calculated value with the prescribed quadrature traction before the
surface traction enters the momentum residual.

The implemented k-exact treatment is consequently:

- `faceGradStencil()` is empty for a `solidTraction` face;
- `faceGradCoeffs()` is empty for that face;
- `fGrad()` leaves its temporary face gradient equal to zero;
- `hofvm` does not add reconstruction columns for that face.

This zero is not a zero-traction boundary condition. It is an unused temporary
gradient that is overwritten by the prescribed traction.

#### Prescribed-displacement boundary reconstruction

This section first gives the idea in scalar notation. For displacement, the
same geometry coefficients are applied separately to its x, y and z
components.

A `fixedDisplacement` value is a value at the boundary face. The unknown in
`kExactLeastSquares` is a cell average. These two values have different
meanings, so the prescribed face value cannot be inserted into the stencil as
if it were another cell average. Instead, it is added as a weighted point-value
observation in the least-squares reconstruction of the adjacent cell.

##### Notation

For one boundary cell, let:

- \(P\) denote the owner cell;
- \(N\) denote one cell in the normal cell stencil \(\mathcal{S}_P\);
- \(r\) denote one prescribed boundary quadrature point attached to \(P\);
- \(q\) denote the quadrature point where the gradient is required;
- \(\overline{u}_P\) and \(\overline{u}_N\) denote cell averages;
- \(u_{D,r}\) denote the prescribed value at boundary point \(r\);
- \(\boldsymbol{a}_P\) denote the scaled polynomial derivatives in cell \(P\).

The letters \(r\) and \(q\) can identify the same point, but do not have to.
The gradient at one point \(q\) can depend on all prescribed points \(r\)
attached to the same owner cell.

##### Step 1: retain the mean-free cell polynomial

The boundary reconstruction uses the same mean-free basis as the interior
cell reconstruction:

$$
u_P^R(\boldsymbol{x})
=
\overline{u}_P
+
\sum_m \psi_{P,m}(\boldsymbol{x})a_{P,m},
$$

where one basis function associated with exponent
\(\boldsymbol{\alpha}_m=(i,j,k)\) is

$$
\psi_{P,m}(\boldsymbol{x})
=
\frac{
(\boldsymbol{x}-\boldsymbol{x}_P)^{\boldsymbol{\alpha}_m}
-M_P^{\boldsymbol{\alpha}_m}
}
{i!j!k!\,h_P^{i+j+k}}.
$$

Its volume average over \(P\) is zero. Therefore,

$$
\frac{1}{V_P}\int_{V_P}u_P^R\,\mathrm{d}V
=\overline{u}_P
$$

without adding a separate cell-average constraint. The constant polynomial
coefficient is still eliminated exactly as in `calcCellCoeffs()`.

##### Step 2: combine neighbour and boundary equations

The normal neighbour equations are

$$
\boldsymbol{Q}_P^T\boldsymbol{a}_P=\boldsymbol{d}_P,
\qquad
d_{P,N}=\overline{u}_N-\overline{u}_P.
$$

These equations are fitted in the weighted least-squares sense. Each
fixed-displacement quadrature point supplies a different kind of equation:

$$
\boldsymbol{D}_P\boldsymbol{a}_P=\boldsymbol{b}_P,
\qquad
b_{P,r}=u_{D,r}-\overline{u}_P.
$$

Row \(r\) of \(\boldsymbol{D}_P\) is simply the mean-free polynomial basis
evaluated at the boundary point \(\boldsymbol{x}_r\):

$$
(D_P)_{r,m}=\psi_{P,m}(\boldsymbol{x}_r).
$$

For example, the two-dimensional quadratic row is

$$
\begin{aligned}
\boldsymbol{D}_{P,r}=
\bigg[&
\frac{r_x}{h_P},
\frac{r_y}{h_P},
\frac{r_x^2-M_{xx}^P}{2h_P^2},
\\
&
\frac{r_xr_y-M_{xy}^P}{h_P^2},
\frac{r_y^2-M_{yy}^P}{2h_P^2}
\bigg],
\end{aligned}
$$

where \(\boldsymbol{r}=\boldsymbol{x}_r-\boldsymbol{x}_P\). This is a point
evaluation of the owner polynomial. It is not a cell-average equation.

The complete mathematical problem is

$$
\boxed{
\underset{\boldsymbol{a}_P}{\operatorname{minimise}}\quad
(\boldsymbol{Q}_P^T\boldsymbol{a}_P-\boldsymbol{d}_P)^T
\boldsymbol{W}_P
(\boldsymbol{Q}_P^T\boldsymbol{a}_P-\boldsymbol{d}_P)
+
(\boldsymbol{D}_P\boldsymbol{a}_P-\boldsymbol{b}_P)^T
\boldsymbol{W}_P^D
(\boldsymbol{D}_P\boldsymbol{a}_P-\boldsymbol{b}_P).
}
$$

The diagonal matrix \(\boldsymbol{W}_P\) contains the same neighbour weights
used by `calcCellCoeffs()`. The diagonal matrix
\(\boldsymbol{W}_P^D\) contains the spatial weights of the boundary
quadrature points. Thus, boundary values are observations in the same weighted
fit as the neighbouring cell averages; they are not exact equality constraints
and no arbitrarily large penalty is used.

##### Step 3: collect every boundary point belonging to the owner cell

`calcFaceGradCoeffs()` first loops over all patches selected by
`includePatchInStencils_`. For a displacement reconstruction, this mask is
formed from `D.boundaryField()[patchI].fixesValue()`.

For every selected boundary face, every face quadrature point is appended to
the list of its owner cell. One address contains

```text
{meshFaceID, quadraturePointID}
```

Here `meshFaceID` is the face index in the local mesh, not the face index
within its boundary patch.

Suppose cell \(P\) touches two fixed-displacement faces. The boundary-data list
for \(P\) contains the points from both faces. When the gradient is evaluated
on either face, one cell polynomial is fitted to all those boundary values. A
boundary observation is never discarded merely because it belongs to the
other face.

##### Step 4: start from the unconstrained reconstruction map

The unconstrained weighted least-squares Hessian and reconstruction map are

$$
\boldsymbol{H}_P
=
\boldsymbol{Q}_P\boldsymbol{W}_P\boldsymbol{Q}_P^T,
$$

$$
\boldsymbol{A}_P
=
\boldsymbol{H}_P^{-1}\boldsymbol{Q}_P\boldsymbol{W}_P.
$$

Therefore, without boundary data,

$$
\boldsymbol{a}_P^{\mathrm{free}}
=
\boldsymbol{A}_P\boldsymbol{d}_P.
$$

`calcCellCoeffs()` already calculated \(\boldsymbol{A}_P\). It stored its rows
as gradient, second-derivative and third-derivative coefficients after removing
the factors \(h_P\), \(h_P^2\) and \(h_P^3\). The boundary code reconstructs
\(\boldsymbol{A}_P\) by collecting those rows and restoring the corresponding
power of \(h_P\). It does not solve the complete neighbour system again.

The code also needs \(\boldsymbol{H}_P^{-1}\). It obtains it from the identity

$$
\boxed{
\boldsymbol{H}_P^{-1}
=
\boldsymbol{A}_P\boldsymbol{W}_P^{-1}\boldsymbol{A}_P^T.
}
$$

To see why this works, substitute the definition of \(\boldsymbol{A}_P\):

$$
\begin{aligned}
\boldsymbol{A}_P\boldsymbol{W}_P^{-1}\boldsymbol{A}_P^T
&=
\boldsymbol{H}_P^{-1}
\boldsymbol{Q}_P\boldsymbol{W}_P
\boldsymbol{W}_P^{-1}
\boldsymbol{W}_P\boldsymbol{Q}_P^T
\boldsymbol{H}_P^{-1}
\\
&=
\boldsymbol{H}_P^{-1}
\boldsymbol{H}_P
\boldsymbol{H}_P^{-1}
\\
&=
\boldsymbol{H}_P^{-1}.
\end{aligned}
$$

##### Step 5: apply the weighted boundary update

Each boundary point uses the same spatial weight function as a neighbouring
cell. For boundary point \(r\),

$$
d_{P,r}^{D}=|\boldsymbol{x}_r-\boldsymbol{x}_P|,
\qquad
w_{P,r}^{D}
=
\texttt{weightFunc().weight}
\left(d_{P,r}^{D},d_P^{\max}\right),
$$

where \(d_P^{\max}\) is the largest cell-centre distance in the normal cell
stencil. These weights form

$$
\boldsymbol{W}_P^D
=
\operatorname{diag}
\left(w_{P,1}^{D},\ldots,w_{P,N_c}^{D}\right).
$$

This explains the difference from `movingLeastSquares`. That reconstruction
is centred directly at a boundary quadrature point, so the distance of that
point from the reconstruction origin is zero and its spatial weight is one.
The k-exact polynomial is centred at the owner-cell centre, so the actual
non-zero distance from the cell centre to the boundary point is used.

Adding the boundary rows to the least-squares problem gives

$$
\widetilde{\boldsymbol{H}}_P
=
\boldsymbol{H}_P
+
\boldsymbol{D}_P^T\boldsymbol{W}_P^D\boldsymbol{D}_P,
$$

$$
\widetilde{\boldsymbol{g}}_P
=
\boldsymbol{Q}_P\boldsymbol{W}_P\boldsymbol{d}_P
+
\boldsymbol{D}_P^T\boldsymbol{W}_P^D\boldsymbol{b}_P,
$$

and

$$
\widetilde{\boldsymbol{H}}_P\boldsymbol{a}_P
=
\widetilde{\boldsymbol{g}}_P.
$$

The implementation does not rebuild and factorise the full augmented system.
It reuses the already calculated cell-only map
\(\boldsymbol{A}_P\) and inverse Hessian
\(\boldsymbol{H}_P^{-1}\). Define

$$
\boldsymbol{S}_P
=
(\boldsymbol{W}_P^D)^{-1}
+
\boldsymbol{D}_P\boldsymbol{H}_P^{-1}\boldsymbol{D}_P^T,
$$

and

$$
\boxed{
\boldsymbol{K}_P
=
\boldsymbol{H}_P^{-1}\boldsymbol{D}_P^T\boldsymbol{S}_P^{-1}.
}
$$

`calcFaceGradCoeffs()` solves systems containing \(\boldsymbol{S}_P\) with
Eigen's column-pivoted Householder QR decomposition. The Woodbury identity then
gives the augmented weighted least-squares solution as

$$
\boxed{
\boldsymbol{a}_P
=
\boldsymbol{a}_P^{\mathrm{free}}
+
\boldsymbol{K}_P
(\boldsymbol{b}_P
-\boldsymbol{D}_P\boldsymbol{a}_P^{\mathrm{free}}).
}
$$

The term in parentheses is the boundary-value residual of the cell-only
polynomial. The matrix \(\boldsymbol{K}_P\) balances correction of this
residual against the fit to the neighbouring cell averages. Therefore, the
resulting polynomial generally does not pass through every prescribed point
exactly; it minimises the combined weighted error.

After inserting
\(\boldsymbol{a}_P^{\mathrm{free}}=\boldsymbol{A}_P\boldsymbol{d}_P\),
the result becomes

$$
\boxed{
\boldsymbol{a}_P
=
\underbrace{
(\boldsymbol{I}-\boldsymbol{K}_P\boldsymbol{D}_P)
\boldsymbol{A}_P
}_{\boldsymbol{R}_P}
\boldsymbol{d}_P
+
\underbrace{\boldsymbol{K}_P}_{\boldsymbol{R}_P^D}
\boldsymbol{b}_P.
}
$$

The code names these two maps `cellMap` and `boundaryMap`:

$$
\texttt{cellMap}=\boldsymbol{R}_P,
\qquad
\texttt{boundaryMap}=\boldsymbol{R}_P^D.
$$

The positive diagonal term \((\boldsymbol{W}_P^D)^{-1}\) makes
\(\boldsymbol{S}_P\) invertible even when some boundary rows are linearly
dependent or when there are more boundary points than polynomial derivative
unknowns. Consequently, corner cells with many face quadrature points no
longer fail because of an over-constrained exact system.

The weighted treatment still reproduces a polynomial of degree no greater
than \(p\): for such a polynomial, both the neighbouring-cell residuals and
the boundary-point residuals are zero for its exact derivative vector.

##### Step 6: differentiate the weighted polynomial at point q

For every target quadrature point \(q\), the code constructs a matrix
\(\boldsymbol{L}_{P,q}\) that differentiates each basis function at
\(\boldsymbol{x}_q\):

$$
\nabla u_P^R(\boldsymbol{x}_q)
=
\boldsymbol{L}_{P,q}\boldsymbol{a}_P.
$$

For the two-dimensional quadratic polynomial, this operation is

$$
\begin{aligned}
\frac{\partial u}{\partial x}(\boldsymbol{x}_q)
&=
\frac{a_x}{h_P}
+\frac{a_{xx}r_{q,x}}{h_P^2}
+\frac{a_{xy}r_{q,y}}{h_P^2},
\\
\frac{\partial u}{\partial y}(\boldsymbol{x}_q)
&=
\frac{a_y}{h_P}
+\frac{a_{xy}r_{q,x}}{h_P^2}
+\frac{a_{yy}r_{q,y}}{h_P^2},
\end{aligned}
$$

where \(\boldsymbol{r}_q=\boldsymbol{x}_q-\boldsymbol{x}_P\). Cubic terms add
the expected quadratic dependence on \(\boldsymbol{r}_q\), and the
three-dimensional version also contains z derivatives.

Multiplication by the two reconstruction maps gives

$$
\boldsymbol{C}_{P,q}^{\mathrm{cell}}
=
\boldsymbol{L}_{P,q}\boldsymbol{R}_P,
\qquad
\boldsymbol{C}_{P,q}^{D}
=
\boldsymbol{L}_{P,q}\boldsymbol{R}_P^D.
$$

The first matrix contains one vector coefficient per stencil cell. The second
contains one vector coefficient per prescribed boundary point.

##### Step 7: convert differences into absolute coefficients

Before storage, the gradient has the difference form

$$
\nabla u_f(\boldsymbol{x}_q)
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{N}^{q}
(\overline{u}_N-\overline{u}_P)
+
\sum_{r\in\mathcal{Q}_P^D}
\boldsymbol{C}_{D,r}^{q}
(u_{D,r}-\overline{u}_P).
$$

`fGrad()` is easier to evaluate using direct multiplication by absolute cell
and boundary values. Expanding the differences gives

$$
\begin{aligned}
\nabla u_f(\boldsymbol{x}_q)
={}&
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{N}^{q}\overline{u}_N
\\
&-
\left(
\sum_{N\in\mathcal{S}_P}\boldsymbol{C}_{N}^{q}
+
\sum_{r\in\mathcal{Q}_P^D}\boldsymbol{C}_{D,r}^{q}
\right)\overline{u}_P
\\
&+
\sum_{r\in\mathcal{Q}_P^D}
\boldsymbol{C}_{D,r}^{q}u_{D,r}.
\end{aligned}
$$

Therefore, the owner-cell coefficient stored in `faceGradCoeffs()` is

$$
\boxed{
\boldsymbol{C}_{P}^{q}
=
-\sum_{N\in\mathcal{S}_P}\boldsymbol{C}_{N}^{q}
-\sum_{r\in\mathcal{Q}_P^D}\boldsymbol{C}_{D,r}^{q}.
}
$$

This is why `calcFaceGradCoeffs()` adds the neighbour and boundary
coefficients to `coefficientSum` and subtracts that sum from the owner
position.

If every cell average and every prescribed value is the same constant, all
terms cancel and the boundary gradient is exactly zero. This constant-field
property is tested explicitly.

##### Step 8: store unknown and prescribed contributions separately

The implemented storage is:

- `faceGradStencilPtr_` stores the global IDs of the owner cell and all cells
  in its reconstruction stencil. The public `faceGradStencil()` accessor
  exposes this list to `fGrad()` and `hofvm`.
- `faceGradCoeffsPtr_` stores the vector coefficient multiplying every cell
  average at each target quadrature point. The public `faceGradCoeffs()`
  accessor exposes these coefficients.
- `faceBoundaryDataStencilPtr_` privately stores every
  `{meshFaceID, quadraturePointID}` supplying prescribed data to that owner
  reconstruction.
- `faceBoundaryDataCoeffsPtr_` privately stores the vector coefficient
  multiplying every prescribed value at every target point.

The prescribed boundary values deliberately do not appear in
`faceGradStencil()`. They are known data, not global cell unknowns.

##### Step 9: evaluate fGrad and assemble the Jacobian

For a vector boundary of type `fixedDisplacement` or `fixedRotation`,
`patchFaceQuadValues()` calls its virtual `evaluateQuadrature()` function.
The base `fixedDisplacement` implementation repeats the current face value at
all quadrature points of that face. `fixedRotation` evaluates the rigid-rotation
displacement separately at every quadrature point. A derived
`fixedDisplacement` condition can also override the function and return
distinct values; the manufactured-solution boundary condition uses this route.

The templated `fGrad()` then performs two sums:

1. public face coefficients multiplied by local or remotely collected cell
   averages;
2. private boundary coefficients multiplied by prescribed boundary values.

The explicit high-order residual therefore contains the complete weighted
boundary gradient. `hofvm` sees only `faceGradStencil()` and
`faceGradCoeffs()`, so its
Jacobian contains derivatives with respect to cell unknowns only. It does not
create matrix columns for prescribed values. Their effect remains on the
residual side of the nonlinear equation.

The weighted-boundary validation reported for this change uses
`highOrderJacobian false`. A focused repeat with `highOrderJacobian true`
remains useful when the high-order Jacobian work is resumed.

##### Boundary equations mapped to the code

| Operation | Implementation |
|---|---|
| Mark prescribed-value patches | `fixedValuePatchMask()` in `solidModel.C` |
| Store the mask in the scheme | `includePatchInStencils_` |
| Build owner plus neighbour addressing | `makeFaceGradStencil()` |
| Collect cell boundary points | `cellBoundaryData` in `calcFaceGradCoeffs()` |
| Reconstruct \(\boldsymbol{A}_P\) | `scaledDerivativeCoeff` lambda |
| Evaluate owner moments in \(\boldsymbol{D}_P\) | `ownerCentralMoment` lambda |
| Evaluate boundary spatial weights | `boundaryWinv` |
| Form \(\boldsymbol{S}_P\) and solve it | `boundaryQr` |
| Store \(\boldsymbol{R}_P\) | `cellMap` |
| Store \(\boldsymbol{R}_P^D\) | `boundaryMap` |
| Differentiate at point \(q\) | matrix `L` |
| Store unknown-cell gradient terms | `faceGradCoeffsPtr_` |
| Store prescribed-value gradient terms | `faceBoundaryDataCoeffsPtr_` |
| Read current prescribed values | `patchFaceQuadValues()` |
| Evaluate both coefficient sums | templated `fGrad()` |
| Assemble unknown-cell Jacobian columns | `hofvm` |

##### Lazy construction and mesh changes

Boundary addressing and coefficients are created lazily together with the
other face coefficients. `clear()` releases
`faceBoundaryDataStencilPtr_` and `faceBoundaryDataCoeffsPtr_` along with the
cell, face, stencil and quadrature caches. They are therefore rebuilt from the
new geometry when the reconstruction is next requested after a mesh update.

##### Current limitations

The present boundary implementation has deliberate limits:

- weighted boundary reconstruction is implemented for vector
  `fixedDisplacement` and `fixedRotation` conditions;
- the scalar fixed-value `patchFaceQuadValues()` overload is not implemented;
- another vector boundary condition for which `fixesValue()` is true causes a
  clear run-time error;
- the base `fixedDisplacement` condition repeats one face value at every
  quadrature point; spatially varying conditions must override
  `evaluateQuadrature()`;
- linearly dependent or overdetermined boundary rows are accepted because
  they are weighted observations rather than exact constraints;
- symmetry reconstruction is still to be implemented;
- the focused fixed-displacement polynomial test currently runs in serial;
  processor-face reconstruction is tested separately in parallel.

The fixed-displacement data modifies only the boundary-face reconstruction.
It does not modify the permanent cell-centred derivative coefficients and does
not change the cell-average unknown used by the time scheme.

#### Cell-coefficient equations mapped to the code

| Mathematical operation | C++ implementation |
|---|---|
| Generate \(x,y,x^2,xy,\ldots\) | **generateExponents()** |
| Locate derivative rows | **rowOf()**, **calcDerivativeRows()** |
| Read \(M_{xx},M_{xy},M_{xxy},\ldots\) | **centralMoment()** |
| Expand an owner-centred term over \(N\) | **averageMonomial()** |
| Subtract the owner moment | **neighbourAverage - ownerAverage** |
| Apply factorial and length scaling | assignment to **Q(p, cI)** |
| Compute the geometry-to-derivative map | **QRSolve(Q, W)** |
| Remove derivative scaling | rows of **A** divided by \(h,h^2,h^3\) |
| Apply map to averages | **grad()**, **secondGrad()**, **thirdGrad()** |

When checking a new term such as **xxy**, trace it through these nine
operations. This separates the moment expansion, factorial, length scaling,
row selection, and field evaluation.

### Compact notation used by the generic implementation

For an owner cell (P), define

$$
\boldsymbol{r}_P = \boldsymbol{x}-\boldsymbol{x}_P,
\qquad
M_P^{\boldsymbol{\alpha}}
=
\frac{1}{V_P}
\int_{V_P}\boldsymbol{r}_P^{\boldsymbol{\alpha}}\,\mathrm{d}V.
$$

The implementation uses the factorial-scaled Taylor form

$$
\phi_P^R(\boldsymbol{x})
=
\overline{\phi}_P
+
\sum_{1\le |\boldsymbol{\alpha}|\le p}
\frac{D^{\boldsymbol{\alpha}}\phi_P}
     {\boldsymbol{\alpha}!}
\left(
\boldsymbol{r}_P^{\boldsymbol{\alpha}}
-M_P^{\boldsymbol{\alpha}}
\right).
$$

Subtracting the owner moment from each basis function makes the volume average
of the reconstruction over cell \(P\) exactly equal to
\(\overline{\phi}_P\). Consequently, the constant polynomial coefficient does
not appear among the least-squares unknowns.

The exponent list generated by `generateExponents()` contains only the
non-constant terms. This differs intentionally from `movingLeastSquares`,
where the constant coefficient is an unknown and the exponent list starts
with `(0,0,0)`. In k-exact reconstruction the owner-cell average already fixes
the constant contribution, so no constant row is assembled or solved.

The implementation assumes that first central moments are zero because the
OpenFOAM cell centre is the volume centroid:

$$
M_P^{100}=M_P^{010}=M_P^{001}=0.
$$

This assumption should be checked numerically using the same volume quadrature
employed for the higher moments.

### Characteristic-length scaling

For each owner-cell stencil, the code defines

$$
h_P = 2\max_{N\in\mathcal{S}_P}
\left|\boldsymbol{x}_N-\boldsymbol{x}_P\right|.
$$

The unknown solved by the scaled least-squares system is

$$
a_P^{\boldsymbol{\alpha}}
=h_P^{|\boldsymbol{\alpha}|}
D^{\boldsymbol{\alpha}}\phi_P.
$$

After solving, `calcCellCoeffs()` divides first-, second- and third-derivative
rows by (h_P), (h_P^2) and (h_P^3), respectively.

### Neighbour cell-average equations

Let

$$
\boldsymbol{d}_{PN}=\boldsymbol{x}_N-\boldsymbol{x}_P.
$$

Inside neighbour cell (N),

$$
\boldsymbol{r}_P=\boldsymbol{d}_{PN}+\boldsymbol{r}_N.
$$

The neighbour average of an owner-centred monomial is evaluated with the
multi-index binomial expansion

$$
\left\langle
\boldsymbol{r}_P^{\boldsymbol{\alpha}}
\right\rangle_N
=
\sum_{\boldsymbol{\beta}\le\boldsymbol{\alpha}}
{\boldsymbol{\alpha}\choose\boldsymbol{\beta}}
\boldsymbol{d}_{PN}^{\boldsymbol{\alpha}-\boldsymbol{\beta}}
M_N^{\boldsymbol{\beta}}.
$$

#### Complete quadratic expansion of averageMonomial()

The general three-loop expression is derived above. For quadratic monomials,
the individual loop combinations can be written explicitly. This makes clear
which offset and central-moment term is generated by each loop.

If an exponent is zero, its loop executes only once. For example, the `xx`
term has \((i,j,k)=(2,0,0)\). Therefore, the `b` and `c` loops each have only
the value zero, while the `a` loop generates the three terms

$$
\begin{aligned}
\left\langle r_{P,x}^{2}\right\rangle_N
={}&
{2\choose0}d_x^2M_N^{000}
+{2\choose1}d_xM_N^{100}
+{2\choose2}M_N^{200}
\\
={}&d_x^2+2d_xM_x^N+M_{xx}^N.
\end{aligned}
$$

For the mixed `xy` term, \((i,j,k)=(1,1,0)\). The `a` and `b` loops generate
four combinations, while the `c` loop again executes once:

$$
\begin{aligned}
\left\langle r_{P,x}r_{P,y}\right\rangle_N
={}&
d_xd_yM_N^{000}
+d_yM_N^{100}
+d_xM_N^{010}
+M_N^{110}
\\
={}&d_xd_y+d_yM_x^N+d_xM_y^N+M_{xy}^N.
\end{aligned}
$$

The complete set of three-dimensional quadratic averages before setting the
first central moments to zero is

$$
\begin{aligned}
\left\langle r_{P,x}^{2}\right\rangle_N
&=d_x^2+2d_xM_x^N+M_{xx}^N,\\
\left\langle r_{P,x}r_{P,y}\right\rangle_N
&=d_xd_y+d_yM_x^N+d_xM_y^N+M_{xy}^N,\\
\left\langle r_{P,x}r_{P,z}\right\rangle_N
&=d_xd_z+d_zM_x^N+d_xM_z^N+M_{xz}^N,\\
\left\langle r_{P,y}^{2}\right\rangle_N
&=d_y^2+2d_yM_y^N+M_{yy}^N,\\
\left\langle r_{P,y}r_{P,z}\right\rangle_N
&=d_yd_z+d_zM_y^N+d_yM_z^N+M_{yz}^N,\\
\left\langle r_{P,z}^{2}\right\rangle_N
&=d_z^2+2d_zM_z^N+M_{zz}^N.
\end{aligned}
$$

Because the neighbour cell centre is its volume centroid,

$$
M_x^N=M_y^N=M_z^N=0,
$$

so the expressions actually returned for the quadratic basis are

$$
\begin{aligned}
\left\langle r_{P,x}^{2}\right\rangle_N &= d_x^2+M_{xx}^N,\\
\left\langle r_{P,x}r_{P,y}\right\rangle_N &= d_xd_y+M_{xy}^N,\\
\left\langle r_{P,x}r_{P,z}\right\rangle_N &= d_xd_z+M_{xz}^N,\\
\left\langle r_{P,y}^{2}\right\rangle_N &= d_y^2+M_{yy}^N,\\
\left\langle r_{P,y}r_{P,z}\right\rangle_N &= d_yd_z+M_{yz}^N,\\
\left\langle r_{P,z}^{2}\right\rangle_N &= d_z^2+M_{zz}^N.
\end{aligned}
$$

After subtracting the corresponding owner moment and applying the
factorial-based Taylor scaling, the six quadratic matrix entries are

$$
\begin{aligned}
Q_{xx,N} &=
\frac{d_x^2+M_{xx}^N-M_{xx}^P}{2h_P^2},\\
Q_{xy,N} &=
\frac{d_xd_y+M_{xy}^N-M_{xy}^P}{h_P^2},\\
Q_{xz,N} &=
\frac{d_xd_z+M_{xz}^N-M_{xz}^P}{h_P^2},\\
Q_{yy,N} &=
\frac{d_y^2+M_{yy}^N-M_{yy}^P}{2h_P^2},\\
Q_{yz,N} &=
\frac{d_yd_z+M_{yz}^N-M_{yz}^P}{h_P^2},\\
Q_{zz,N} &=
\frac{d_z^2+M_{zz}^N-M_{zz}^P}{2h_P^2}.
\end{aligned}
$$

The diagonal entries contain the factor \(2!=2\). The mixed entries contain
\(1!1!=1\). In two dimensions, all expressions containing z are omitted, so
only `xx`, `xy`, and `yy` remain.

`centralMoment()` supplies \(M_N^{\boldsymbol{\beta}}\): order zero is one,
order one is assumed zero, and orders two and three come from
`fvMeshQuadrature`. `averageMonomial()` performs the binomial sum. This generic
calculation is used for all supported exponents instead of writing separate
expressions for `xx`, `xy`, `xxy`, and so on.

For every non-constant exponent and stencil neighbour, the matrix entry is

$$
Q_{\boldsymbol{\alpha}N}
=
\frac{
\left\langle
\boldsymbol{r}_P^{\boldsymbol{\alpha}}
\right\rangle_N
-M_P^{\boldsymbol{\alpha}}
}
{h_P^{|\boldsymbol{\alpha}|}\boldsymbol{\alpha}!}.
$$

For example,

$$
Q_{xx,N}
=
\frac{d_x^2+M_{xx}^N-M_{xx}^P}{2h_P^2},
$$

and

$$
Q_{xxx,N}
=
\frac{
d_x^3+3d_xM_{xx}^N+M_{xxx}^N-M_{xxx}^P
}{6h_P^3}.
$$

The factors two and six arise from the Taylor factorials. If a derivation uses
raw polynomial coefficients such as \(c_{xx}r_x^2\), then
\(D_{xx}=2c_{xx}\). Similarly, \(D_{xxx}=6c_{xxx}\). Mixed-term factors follow
the multi-index factorial; for example, `xy` has factor one and `xxy` has
factor two.

Collecting all neighbours gives

$$
\boldsymbol{Q}_P^T\boldsymbol{a}_P=\boldsymbol{b}_P,
\qquad
(b_P)_N=\overline{\phi}_N-\overline{\phi}_P.
$$

In the code, `Q` is stored transposed relative to the usual equation layout:
its rows are non-constant basis terms and its columns are stencil neighbours.
The numbers of solved basis terms are:

| Order | Two dimensions | Three dimensions |
|------:|---------------:|-----------------:|
| 1 | 2 | 3 |
| 2 | 5 | 9 |
| 3 | 9 | 19 |

`minNn()` returns exactly these numbers because only the non-constant
coefficients are unknown. `leastSquaresStencil` expects its requested size to
include the central cell, so `makeStencils()` passes
`minNn() + cellStencilExtraCells + 1`. Since `cellsStencil()` excludes the
central cell, the resulting least-squares system has
`minNn() + cellStencilExtraCells` neighbour equations.

### Weighted QR solution

For neighbour weights (w_N), `QRSolve()` forms

$$
\widehat{\boldsymbol{Q}}
=\boldsymbol{Q}\boldsymbol{W}^{1/2},
\qquad
\boldsymbol{W}=\operatorname{diag}(w_N),
$$

and uses a Householder QR decomposition. This corresponds to minimizing

$$
\sum_{N\in\mathcal{S}_P}
w_N
\left(
(\boldsymbol{Q}_P^T\boldsymbol{a}_P)_N-(b_P)_N
\right)^2.
$$

The resulting matrix `A` has one row per non-constant basis term and one column
per stencil neighbour, and satisfies

$$
\boldsymbol{a}_P=\boldsymbol{A}_P\boldsymbol{b}_P.
$$

The same optional condition-number calculation used by
`movingLeastSquares` is retained. Its output field is named
`kExactCellConditionNumber`.

### Stored coefficients and evaluation

Each coefficient list has exactly the same length as the corresponding cell
stencil. There is no additional owner-cell entry. For example, the stored
gradient coefficient is

$$
\boldsymbol{C}_{PN}^{\nabla}
=
\left(
\frac{A_{xN}}{h_P},
\frac{A_{yN}}{h_P},
\frac{A_{zN}}{h_P}
\right),
$$

and `grad()` evaluates

$$
\nabla\phi_P
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{PN}^{\nabla}
\left(
\overline{\phi}_N-\overline{\phi}_P
\right).
$$

`secondGrad()` and `thirdGrad()` use the same neighbour-minus-owner form with
rows of `A` divided by (h_P^2) and (h_P^3). Scalar and vector `grad()`
share a templated implementation. In two dimensions, all derivative tensor
components containing z are explicitly set to zero.

### `valueAtPoint`

The common reconstruction interface evaluates a scalar or vector field from a
specified cell at an arbitrary point:

```cpp
reconstruction.valueAtPoint(field, cellID, x);
```

For `movingLeastSquares`, the stored value is interpreted as the value at the
cell centre. With

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P,
$$

the scalar reconstruction is

$$
u^R(\boldsymbol{x})
=u_P
+\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}\nabla u_P
+\frac{1}{2}\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{H}_P\boldsymbol{r}
+\frac{1}{6}\boldsymbol{T}_P:\boldsymbol{r}^3.
$$

Here, \(\boldsymbol{H}_P\) and \(\boldsymbol{T}_P\) are the reconstructed
second- and third-derivative tensors. Terms above the selected polynomial order
are omitted.

For `kExactLeastSquares`, the stored value is the cell average. The mean-free
basis therefore gives

$$
\begin{aligned}
u^R(\boldsymbol{x})={}&\overline{u}_P
+\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}\nabla u_P\\
&+\frac{1}{2}
\left(
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}\boldsymbol{H}_P\boldsymbol{r}
-\boldsymbol{H}_P:\boldsymbol{M}^{(2)}_P
\right)\\
&+\frac{1}{6}
\left(
\boldsymbol{T}_P:\boldsymbol{r}^3
-\boldsymbol{T}_P:\boldsymbol{M}^{(3)}_P
\right).
\end{aligned}
$$

The moment corrections are essential: without them, the expansion would treat
the cell average as though it were a point value at the centroid. For a vector
field, the same equation is applied independently to each component.

`valueAtPoint()` does not construct `grad()`, `secondGrad()` or `thirdGrad()`
fields over the complete mesh. Instead, it takes only the coefficient row for
the selected cell, contracts those coefficients with the displacement from the
cell centre to the requested point, and multiplies the resulting scalar
coefficients by the values in that cell's stencil. The same calculation works
for scalar and vector fields.

The coefficient tables themselves are still mesh-wide cached data. They are
normally already available to the high-order solver and are only constructed
on their first use or after a mesh change.

In a parallel run, the operation must be called by every processor because the
selected cell stencil may contain remote cells. The processor that owns the
requested cell passes its local `cellID`; all other processors pass `-1` and
receive zero. All processors participate in the remote-field exchange, but
only the owning processor evaluates a cell reconstruction. The resulting value
can then be combined with a global sum. `solidPointDisplacement` follows this
pattern.

For a parallel stencil, neighbour cell averages, cell centres, and required
second- or third-order moments are obtained using the communication functions
of `leastSquaresStencil`. Second moments are requested only for (p\ge2), and
third moments only for (p\ge3).

#### `cellValueCoeffsAtPoint`

The internal virtual helper has the following form:

```cpp
void cellValueCoeffsAtPoint
(
    const label cellID,
    const point& x,
    UList<scalar>& coeffs
) const;
```

It converts the derivative coefficient rows already stored for cell $P$
into scalar value coefficients at an arbitrary point
$\boldsymbol{x}$. It does not read or evaluate a field. Its result
satisfy
$$
u_P^R(\boldsymbol{x})
=
A_{PP}(\boldsymbol{x})u_P
+
\sum_{j\in\mathcal{S}_P}
A_{Pj}(\boldsymbol{x})u_j.
$$

Here, \(\mathcal{S}_P\) is the cell stencil. The output list has
`stencil().cellsStencil()[cellID].size() + 1` entries. The stencil-cell
coefficients retain the stencil ordering, and the last entry is always
be the owner-cell coefficient:

$$
\left[
A_{P1},A_{P2},\ldots,A_{Pn},A_{PP}
\right].
$$

The coefficients are scalar. Consequently, exactly the same coefficient list
can reconstruct a scalar field or each component of a vector field.

##### Moving-least-squares coefficients

For `movingLeastSquares`, $u_P$ is the field value at the cell centre. Let

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P.
$$

For every coefficient row $j$, including the final central-cell row, define

$$
\begin{aligned}
b_{Pj}(\boldsymbol{x})={}&
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{G}_{Pj}\\
&+\frac{1}{2}
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{H}_{Pj}\boldsymbol{r}\\
&+\frac{1}{6}
\boldsymbol{T}_{Pj}:\boldsymbol{r}^3.
\end{aligned}
$$

In this equation, \(\boldsymbol{G}_{Pj}\),
\(\boldsymbol{H}_{Pj}\), and \(\boldsymbol{T}_{Pj}\) are respectively the
gradient-, second-derivative-, and third-derivative-coefficient rows. Terms
above the selected polynomial order are omitted. The value coefficients are
then

$$
A_{Pj}=b_{Pj},
\qquad j\in\mathcal{S}_P,
$$

and

$$
A_{PP}=1+b_{PP}.
$$

The one in the owner coefficient is the original cell-centre value in the
Taylor expansion. The term $b_{PP}$ is the contribution of the central-cell
derivative coefficient row. Since the derivative coefficient rows annihilate
a constant field, they also satisfy

$$
b_{PP}+\sum_{j\in\mathcal{S}_P}b_{Pj}=0
$$

to numerical tolerance. Consequently, the moving-least-squares value
coefficients also sum to one.

##### K-exact coefficients

For `kExactLeastSquares`, the unknown is the cell average
$\overline{u}_P$. For each stencil cell $j$, define

$$
\begin{aligned}
b_{Pj}(\boldsymbol{x})={}&
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{G}_{Pj}\\
&+\frac{1}{2}
\left(
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{H}_{Pj}\boldsymbol{r}
-\boldsymbol{H}_{Pj}:\boldsymbol{M}^{(2)}_P
\right)\\
&+\frac{1}{6}
\left(
\boldsymbol{T}_{Pj}:\boldsymbol{r}^3
-\boldsymbol{T}_{Pj}:\boldsymbol{M}^{(3)}_P
\right).
\end{aligned}
$$

The owner-cell reconstruction is currently expressed in the
neighbour-minus-owner form
$$
u_P^R(\boldsymbol{x})
=
\overline{u}_P
+
\sum_{j\in\mathcal{S}_P}
b_{Pj}(\boldsymbol{x})
\left(\overline{u}_j-\overline{u}_P\right).
$$

Expanding this expression gives

$$
u_P^R(\boldsymbol{x})
=
\left(1-\sum_{j\in\mathcal{S}_P}b_{Pj}\right)
\overline{u}_P
+
\sum_{j\in\mathcal{S}_P}b_{Pj}\overline{u}_j.
$$

Therefore, `cellValueCoeffsAtPoint()` returns

$$
A_{Pj}=b_{Pj},
\qquad j\in\mathcal{S}_P,
$$

and

$$
A_{PP}=1-\sum_{j\in\mathcal{S}_P}b_{Pj}.
$$

It follows directly that

$$
A_{PP}+\sum_{j\in\mathcal{S}_P}A_{Pj}=1.
$$

Thus, applying the coefficients to a constant field returns exactly that
constant, up to round-off error.

For example, a first-order k-exact reconstruction has

$$
b_{Pj}=\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{G}_{Pj},
$$

and hence

$$
u_P^R(\boldsymbol{x})
=
\left(1-\sum_jb_{Pj}\right)\overline{u}_P
+
\sum_jb_{Pj}\overline{u}_j.
$$

##### How the coefficients are used

Once the required local and remote field values are available, evaluation can
be written schematically as

```cpp
Type value = coeffs.last()*field[cellID];

forAll(cellStencil, stencilI)
{
    value += coeffs[stencilI]*stencilFieldValue[stencilI];
}
```

The helper itself performs only geometry-dependent contractions. It does not
perform processor communication, enforce boundary conditions, integrate
over the cell, or multiply the coefficients by field values. Those operations
belong to the caller.

This makes `cellValueCoeffsAtPoint()` the value counterpart of
`cellGradCoeffsAtPoint()`:

$$
u_P^R(\boldsymbol{x})
=\sum_J A_{PJ}(\boldsymbol{x})u_J,
\qquad
\nabla u_P^R(\boldsymbol{x})
=\sum_J\boldsymbol{C}_{PJ}(\boldsymbol{x})u_J.
$$

The value coefficients are used by `valueAtPoint()` and by the face-centre
evaluation described below. They can also be reused by the future exact
linearisation of the stabilisation term.

#### `faceCentreValues`

`leastSquaresScheme` provides scalar and vector overloads with the common
interface

```cpp
reconstruction.faceCentreValues
(
    field,
    ownerValues,
    neighbourValues
);
```

Both output arguments are surface fields. For an internal face (f), they
contain two independent reconstructions at the same face centre:

$$
u_{L,f}=u_P^R(\boldsymbol{C}_f),
\qquad
u_{R,f}=u_N^R(\boldsymbol{C}_f).
$$

The face-centre coefficient storage belongs to `leastSquaresScheme` because
the addressing and evaluation procedure are common to both reconstruction
types. The scheme-specific `cellValueCoeffsAtPoint()` implementation supplies
each coefficient row.

The coefficients are constructed lazily. On the first call,
`makeFaceCentreValueCoeffs()` creates:

- one owner-side coefficient row for every face;
- one neighbour-side coefficient row for every internal face.

For an internal face, construction calls

```cpp
cellValueCoeffsAtPoint(owner[faceI], Cf[faceI], ownerCoeffs);
cellValueCoeffsAtPoint(neighbour[faceI], Cf[faceI], neighbourCoeffs);
```

For a boundary face, only the local owner-side row exists. All rows use the
same owner-last layout described above. Since the coefficients depend only on
the reconstruction and mesh geometry, the scalar and vector overloads share
the same cache.

Every call to `faceCentreValues()` exchanges the current remote stencil-field
values once and multiplies them by the cached coefficients. The field values
cannot be cached because they change during the solution, but no QR solution,
moment calculation or derivative contraction is repeated.

At a processor face, each processor evaluates its local owner-cell
reconstruction. The completed reconstructed values are then exchanged, giving

$$
u_{L,f}=u_P^R(\boldsymbol{C}_f)
$$

on the local side and

$$
u_{R,f}=u_N^R(\boldsymbol{C}_f)
$$

from the neighbouring processor. A remote processor-local cell label is never
passed to `cellValueCoeffsAtPoint()`.

At an uncoupled physical boundary there is no neighbour-cell reconstruction.
The implementation therefore initializes the neighbour output to the local
owner reconstruction, producing a neutral zero jump. The caller remains
responsible for boundary-condition policy. `alphaStab`, described below,
replaces this neighbour value with the prescribed value on a fixed-value
boundary and retains zero stabilisation on a traction boundary.

Both cached coefficient lists are cleared by
`clearFaceCentreValueCoeffs()`. The derived scheme calls this from its existing
`clear()` function, so mesh motion invalidates the cache and its next use
reconstructs it lazily for the new geometry.

#### `alphaStab`

The high-order `alphaStab` implementation obtains its face-centre values
directly from the selected `leastSquaresScheme`. It no longer constructs a
pointwise Taylor expansion from separately evaluated first, second and third
derivatives. Consequently, the same code is used for `movingLeastSquares` and
`kExactLeastSquares`, while each scheme retains its own value semantics.

For a vector field, `alphaStab` uses `displacementLeastSquares()`. For a scalar
pressure field, it uses `pressureLeastSquares()`. The reconstruction returns
the owner- and neighbour-side values

$$
u_{L,f}=u_P^R(\boldsymbol{C}_f),
\qquad
u_{R,f}=u_N^R(\boldsymbol{C}_f).
$$

On an internal face, the stabilisation is

$$
s_f = \alpha\frac{u_{R,f}-u_{L,f}}
{\max\left(\left|\boldsymbol{n}_f\mathbin{\cdot}
\boldsymbol{d}_f\right|,\epsilon\right)},
$$

where $\boldsymbol{d}_f$ is the existing `deltaVectors` value and $\alpha$ is
the configured `scaleFactor`.

At a processor face, both processors reconstruct their local owner-side value
and `faceCentreValues()` exchanges the completed values. The same jump is then
formed using the vector between the two adjacent cell centres in the
denominator. Derivatives and processor-local cell labels are not exchanged by
`alphaStab`.

On a fixed-value physical boundary, only the owner reconstruction is required.
The right value is the prescribed patch value $u_{b,f}$:

$$
s_f = \alpha\frac{u_{b,f}-u_{L,f}}
{\max\left(\left|\boldsymbol{n}_f\mathbin{\cdot}
(\boldsymbol{C}_f-\boldsymbol{C}_P)\right|,\epsilon\right)}.
$$

This covers `fixedDisplacement` through its `fixesValue()` interface. On
traction, empty, symmetry and other non-fixed physical patches, `alphaStab`
sets the contribution to zero. The implicit `alphaStab` Jacobian is still the
existing Laplacian approximation; exact k-exact stabilisation coefficients are
separate future work.

### Current validation boundary

The current sources and test utilities build with OpenFOAM.com v2412.
`Test-highOrderGrad` reads the reconstruction settings from
`constant/solidProperties` and checks:

- cell-centred gradients of scalar and vector linear, quadratic and cubic
  polynomials;
- scalar and vector internal-face gradients at every face quadrature point;
- scalar and vector gradients at processor-face quadrature points;
- exact constant-field cancellation on all patches whose runtime type is
  exactly `fixedDisplacement`;
- one-sided fixed-displacement gradients for linear, quadratic and cubic
  fields on the first planar patch of that type. Derived analytical boundary
conditions are not replaced by the manufactured test field.

The Cook's membrane regression runs `Test-highOrderGrad` in serial for both
reconstruction types. The spherical-cavity regression runs it with two
processors for both types.

For `kExactLeastSquares`, polynomial cell averages are calculated independently
with third-order cell quadrature. For `movingLeastSquares`, the polynomial is
sampled at the cell centre. In three dimensions, the fields contain every
monomial in x, y and z through the tested order.

The two-dimensional test and the three-dimensional tetrahedral
spherical-cavity test pass through cubic order. In the serial spherical-cavity
test with two fixed-displacement patches, including cells touching both
patches, the maximum constant-field boundary-gradient error is approximately
\(1.1\times10^{-12}\). The maximum fixed-boundary errors for first-, second-
and third-order fields are approximately \(1.1\times10^{-11}\).

The same tetrahedral case passes with a two-processor decomposition. The
latest run exchanged 4601 remote cell values, used remote values in 3872 cell
stencils and 7615 internal-face stencils, and directly tested all 478 processor
faces. Cell, internal-face and processor-face gradient errors through cubic
order were approximately \(10^{-11}\) or smaller. The focused
fixed-displacement test is currently skipped in parallel; fixed physical
boundaries and processor faces are tested independently.

The canonical `fvMeshQuadrature` ordering was additionally checked using four
processors, where every rank had several processor neighbours. The 3D
tetrahedral test passed all 784 processor faces through cubic order, with a
maximum processor-face gradient error of approximately \(1.4\times10^{-11}\).
The 2D plate-hole test passed all 63 processor faces through its selected
quadratic order for both `kExactLeastSquares` and `movingLeastSquares`.

The reconstruction-aware `alphaStab` was checked with the displacement field
on Cook's membrane. Clean serial runs with both `movingLeastSquares` and
`kExactLeastSquares` converged in three nonlinear iterations. The moving
least-squares run reproduced the stored tip-displacement magnitude of
0.0402181. A two-processor, 20-by-20 k-exact Cook case, with 20 processor faces
and both fixed-displacement and traction boundaries, also converged in three
nonlinear iterations.

A temporary second-order Cook's membrane case was also used as a solver-level
integration check. It contains a real `fixedDisplacement` boundary and
converged in three SNES iterations with both `highOrderJacobian false` and
`highOrderJacobian true`. With the analytic high-order Jacobian, the monitored
tip displacements were:

| Reconstruction | \(D_x\) | \(D_y\) | \(|\boldsymbol{D}|\) |
|---|---:|---:|---:|
| `kExactLeastSquares` | -0.0234674 | 0.0317955 | 0.0395180 |
| `movingLeastSquares` | -0.0234752 | 0.0317757 | 0.0395067 |

The displacement-magnitude difference is \(1.13\times10^{-5}\), or about
0.0286 percent relative to the MLS result. This is an integration and
regression check, not a mesh-convergence or reference-solution validation.

The general-polyhedral spherical-cavity mesh still exposes non-zero first
central moments and a mismatch between the sum of cell quadrature weights and
`mesh.V()`. Consequently, its quadrature-generated field values are not exact
cell averages. Remaining work includes:

- resolving general-polyhedral cell-quadrature volume and first-moment errors;
- implementing symmetry reconstruction;
- extending prescribed-value support beyond `fixedDisplacement` and
  `fixedRotation` where needed;
- validating fixed-displacement reconstruction in parallel;
- testing transient and updated-Lagrangian formulations.

## Creating and accessing a reconstruction

Solver-side code first obtains the one manager associated with the mesh:

```cpp
const leastSquaresReconstruction& reconstructions =
    leastSquaresReconstruction::New(mesh);
```

`New(mesh)` searches the mesh object registry. It constructs and registers the
manager on the first call and returns the existing object on later calls. The
mesh owns the manager; the local variable is only a reference.

A solid model can then provide field-specific convenience accessors:

```cpp
const leastSquaresScheme&
solidModel::displacementLeastSquares() const
{
    const leastSquaresReconstruction& reconstructions =
        leastSquaresReconstruction::New(mesh());

    return reconstructions.scheme
    (
        "displacement",
        fixedValuePatchMask(incremental() ? DD_ : D_),
        displacementHighOrderCoeffs()
    );
}
```

The corresponding pressure accessor supplies the pressure dictionary and the
boundary mask obtained from `p()`.

On its first call, `scheme(...)` uses the reconstruction dictionary to create
either `movingLeastSquares` or `kExactLeastSquares`. Later calls return the
stored object. The current manager does not compare a later dictionary or
boundary mask with the data used on the first call. All users of one field role
must therefore request it consistently.

Code using the returned base reference calls common virtual functions:

```cpp
const leastSquaresScheme& reconstruction =
    displacementLeastSquares();

reconstruction.fGrad(D(), gradDQuad);
```

Normal C++ virtual dispatch invokes the implementation belonging to the type
selected from the dictionary.

## Reconstruction dictionary

Displacement and pressure reconstruction must be independently configurable:

```text
highOrderCoeffs
{
    displacement
    {
        type                 movingLeastSquares;
        polynomialOrder      2;
        faceStencilExtraCells 4;
        cellStencilExtraCells 4;

        // Existing MLS settings
        // ...
    }

    pressure
    {
        type                 kExactLeastSquares;
        polynomialOrder      2;
        cellStencilExtraCells 4;

        // k-exact settings
        // ...
    }
}
```

For `kExactLeastSquares`, `cellStencilExtraCells` is preferred. If it is not
specified, the existing `faceStencilExtraCells` entry is used for backward
compatibility.

Omitting `type` selects `movingLeastSquares`, which keeps existing case
dictionaries valid.

The architecture permits a pressure reconstruction to be created, but this
does not by itself implement a high-order mixed displacement-pressure
formulation. Solver residuals, Jacobians and mechanical-law coupling must also
support that formulation before it can be enabled.

## Mesh changes and cached data

Stencil cells, quadrature locations, moments and reconstruction coefficients
depend on mesh geometry. When the mesh moves, the implemented `movePoints()`
callback asks each existing scheme to clear its geometry-dependent data. The
data are then reconstructed lazily on their next access.

For a topology change, field boundary masks may also become invalid. The
implemented `updateMesh()` callback deletes the displacement and pressure
scheme objects. Their boundary masks and geometry data are reconstructed on
the next access.

The migration must include an updated-Lagrangian test because it exercises
this lifecycle; a successful static-mesh case is not sufficient validation of
mesh-object ownership.

## OpenFOAM compatibility

The design uses one type-named mesh object:

```cpp
leastSquaresReconstruction::New(mesh)
```

This form is available in OpenFOAM.com and OpenFOAM.org. The implementation
should not depend on the OpenFOAM.com facility for creating several mesh
objects with user-supplied registry names.

Run-time selection tables, `autoPtr` and the standard type-named mesh-object
lookup provide the common implementation path. Version-specific mesh-motion
callback signatures should follow the compatibility pattern already used by
other solids4foam mesh objects.

## Implementation plan

The refactoring and the new numerical method are being implemented as separate
milestones:

1. Record baseline MLS results for representative serial, parallel and
   moving-mesh cases.
2. Validate the `leastSquaresStencil` rename against the MLS baseline without
   changing its behaviour.
3. Add `leastSquaresScheme`, derive and register `movingLeastSquares`, and
   preserve the existing MLS results.
4. Add the concrete `leastSquaresReconstruction` mesh-object manager and move
   displacement and pressure reconstruction ownership out of `solidModel`.
   This step is complete.
5. Build with supported OpenFOAM.com and OpenFOAM.org versions and repeat the
   serial, parallel and moving-mesh validation. The original framework and
   stencil rename were validated on both versions; the latest boundary changes
   have currently been built only with OpenFOAM.com v2412.
6. Add the registered `kExactLeastSquares` skeleton. This step is complete.
7. Implement its moment data, reconstruction coefficients and
   polynomial-exactness tests. This step is complete for cell, internal-face,
   processor-face and `fixedDisplacement` gradients.
8. Implement symmetry and any additional required prescribed-value boundary
   conditions.
9. Resolve the cell quadrature errors seen on the general-polyhedral mesh.
10. Add transient and updated-Lagrangian tests demonstrating the cell-average
    time-term formulation and mesh-object lifecycle.

Each step should remain numerically equivalent to the previous one until the
new k-exact implementation is deliberately selected.
