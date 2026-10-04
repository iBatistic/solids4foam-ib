/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "quadraticGradientCorrection.H"
#include "kExactLeastSquares.H"
#include "movingLeastSquares.H"
#include <Eigen/Dense>

namespace Foam
{
namespace
{
    // Ordering agrees with symmTensor3rdOrder; mixed terms include factorials.
    const label cubicPowers[10][3] =
    {
        {3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1},
        {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}
    };

    const label secondPowers[6][3] =
    {
        {2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2}
    };
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar quadraticGradientCorrection::monomial
(
    const vector& r,
    const label* powers
)
{
    scalar value = 1;
    for (direction d = 0; d < 3; ++d)
    {
        for (label k = 1; k <= powers[d]; ++k)
        {
            value *= r[d]/k;
        }
    }
    return value;
}


const label* quadraticGradientCorrection::powers(const label term) const
{
    return nTerms_ == 6 ? secondPowers[term] : cubicPowers[term];
}


label quadraticGradientCorrection::cubicTerm
(
    const label i,
    const label j,
    const label k
)
{
    label p[3] = {0, 0, 0};
    ++p[i];
    ++p[j];
    ++p[k];
    for (label term = 0; term < 10; ++term)
    {
        if
        (
            cubicPowers[term][0] == p[0]
         && cubicPowers[term][1] == p[1]
         && cubicPowers[term][2] == p[2]
        )
        {
            return term;
        }
    }
    FatalErrorInFunction << "No cubic term for indices " << i << j << k
        << abort(FatalError);
    return -1;
}


void quadraticGradientCorrection::makeFaceMoments()
{
    const fvMesh& mesh = reconstruction_.mesh();
    const pointField& pts = mesh.points();
    auto& points = compactListListCRef(reconstruction_.quadrature().faceQuadPoints());
    auto& weights = compactListListCRef(reconstruction_.quadrature().faceQuadWeights());

    faceMoment_.setSize(mesh.nFaces(), symmTensor::zero);
    forAll(faceMoment_, faceI)
    {
        const face& f = mesh.faces()[faceI];
        const point& c = mesh.faceCentres()[faceI];

        // Exact second moment about the face centre from the fan of
        // triangles (c, a, b): A*[(ct-c)(ct-c) + sum_v (v-ct)(v-ct)/12]
        symmTensor J(symmTensor::zero);
        scalar area = 0;
        forAll(f, i)
        {
            const point& a = pts[f[i]];
            const point& b = pts[f[f.fcIndex(i)]];
            const scalar At = 0.5*mag((a - c) ^ (b - c));
            const point ct = (a + b + c)/3.0;
            J += At*(sqr(ct - c) + (sqr(a - ct) + sqr(b - ct) + sqr(c - ct))/12.0);
            area += At;
        }

        // Second moment as sampled by the quadrature rule
        forAll(points[faceI], q)
        {
            J -= weights[faceI][q]*sqr(points[faceI][q] - c);
        }

        faceMoment_[faceI] = 0.5*J/area;
    }
}


scalar quadraticGradientCorrection::sample
(
    const label cellI,
    const point& origin,
    const label term
) const
{
    return sampleMonomial(cellI, origin, powers(term));
}


scalar quadraticGradientCorrection::sampleMonomial
(
    const label cellI,
    const point& origin,
    const label* e
) const
{
    const vector r = reconstruction_.mesh().C()[cellI] - origin;
    scalar value = monomial(r, e);

    if (isA<kExactLeastSquares>(reconstruction_))
    {
        const fvMeshQuadrature& quad = momentQuadrature_;
        const vector& m1 = quad.firstOrderCellMoments()[cellI];
        const symmTensor& m2 = quad.secondOrderCellMoments()[cellI];

        for (direction d = 0; d < 3; ++d)
        {
            if (e[d])
            {
                label p[3] = {e[0], e[1], e[2]};
                --p[d];
                value += m1[d]*monomial(r, p);
            }
        }
        for (direction d = 0; d < 6; ++d)
        {
            label p[3];
            bool included = true;
            scalar factorial = 1;
            for (direction axis = 0; axis < 3; ++axis)
            {
                p[axis] = e[axis] - secondPowers[d][axis];
                included = included && p[axis] >= 0;
                if (secondPowers[d][axis] == 2)
                {
                    factorial *= 2;
                }
            }
            if (included)
            {
                value += m2[d]*monomial(r, p)/factorial;
            }
        }
        if (e[0] + e[1] + e[2] == 3)
        {
            scalar factorial = 1;
            for (direction axis = 0; axis < 3; ++axis)
            {
                for (label k = 2; k <= e[axis]; ++k)
                {
                    factorial *= k;
                }
            }
            for (label term = 0; term < 10; ++term)
            {
                if
                (
                    cubicPowers[term][0] == e[0]
                 && cubicPowers[term][1] == e[1]
                 && cubicPowers[term][2] == e[2]
                )
                {
                    value += quad.thirdOrderCellMoments()[cellI][term]/factorial;
                }
            }
        }
    }
    return value;
}


FixedList<scalar, 10> quadraticGradientCorrection::variation
(
    const UList<vector>& firstDerivative,
    const label cellI
) const
{
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());
    auto& weights = compactListListCRef(variationCoeffs_);
    tensor derivative(tensor::zero);
    forAll(stencils[cellI], j)
    {
        derivative += weights[cellI][j]
           *(firstDerivative[stencils[cellI][j]] - firstDerivative[cellI]);
    }
    const symmTensor secondDerivative = symm(derivative);
    FixedList<scalar, 10> result(scalar(0));
    for (direction d = 0; d < 6; ++d)
    {
        result[d] = secondDerivative[d];
    }
    return result;
}


template<class Type>
void quadraticGradientCorrection::makeResponse
(
    const CompactListList<Type>& coefficients
)
{
    const vectorField& C = reconstruction_.mesh().C();
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());
    auto& coeffs = compactListListCRef(coefficients);
    for (label term = 0; term < nTerms_; ++term)
    {
        // At its own expansion centre the exact highest reconstructed
        // derivative of this omitted-degree monomial is zero.
        List<Type> bias(C.size(), pTraits<Type>::zero);
        forAll(stencils, cellI)
        {
            const scalar ownSample = sample(cellI, C[cellI], term);
            forAll(stencils[cellI], j)
            {
                bias[cellI] += coeffs[cellI][j]
                   *(sample(stencils[cellI][j], C[cellI], term) - ownSample);
            }
        }
        forAll(stencils, cellI)
        {
            const FixedList<scalar, 10> response = variation(bias, cellI);
            for (label row = 0; row < nTerms_; ++row)
            {
                responseInverse_[cellI](row,term) =
                    response[row] + (row == term ? 1 : 0);
            }
        }
    }
}


template<class Type>
List<FixedList<scalar, 10>> quadraticGradientCorrection::rawDerivatives
(
    const CompactListList<Type>& coefficients,
    const scalarField& values
) const
{
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());
    auto& coeffs = compactListListCRef(coefficients);
    List<Type> derivative(values.size(), pTraits<Type>::zero);
    forAll(stencils, cellI)
    {
        forAll(stencils[cellI], j)
        {
            derivative[cellI] +=
                coeffs[cellI][j]*(values[stencils[cellI][j]] - values[cellI]);
        }
    }
    List<FixedList<scalar, 10>> result(values.size());
    forAll(result, cellI)
    {
        result[cellI] = variation(derivative, cellI);
    }
    return result;
}


List<FixedList<scalar, 10>> quadraticGradientCorrection::recoveredDerivatives
(
    const scalarField& values
) const
{
    const List<FixedList<scalar, 10>> raw =
        rawDerivatives(reconstruction_.cellGradCoeffs(), values);
    List<FixedList<scalar, 10>> result(values.size());
    forAll(result, cellI)
    {
        result[cellI] = scalar(0);
        for (label i = 0; i < nTerms_; ++i)
        {
            for (label j = 0; j < nTerms_; ++j)
            {
                result[cellI][i] += responseInverse_[cellI](i,j)*raw[cellI][j];
            }
        }
    }
    return result;
}


List<symmTensor> quadraticGradientCorrection::cellSecondDerivatives
(
    const scalarField& values
) const
{
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());
    auto& coeffs = compactListListCRef(reconstruction_.cellSecondGradCoeffs());
    List<symmTensor> secondDerivative(values.size(), symmTensor::zero);
    forAll(stencils, cellI)
    {
        forAll(stencils[cellI], j)
        {
            secondDerivative[cellI] += coeffs[cellI][j]
               *(values[stencils[cellI][j]] - values[cellI]);
        }
    }
    return secondDerivative;
}


List<FixedList<scalar, 10>>
quadraticGradientCorrection::cellThirdDerivatives
(
    const scalarField& values
) const
{
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());
    List<FixedList<scalar, 10>> result(values.size());
    forAll(result, cellI)
    {
        const UList<label>& cells = stencils[cellI];
        const List<FixedList<scalar, 10>>& coeffs = cellCubicCoeffs_[cellI];
        FixedList<scalar, 10>& T = result[cellI];
        T = scalar(0);
        forAll(cells, j)
        {
            for (label m = 0; m < 10; ++m)
            {
                T[m] += coeffs[j][m]*values[cells[j]];
            }
        }
        for (label m = 0; m < 10; ++m)
        {
            T[m] += coeffs[cells.size()][m]*values[cellI];
        }
    }
    return result;
}


List<FixedList<scalar, 10>>
quadraticGradientCorrection::faceThirdDerivatives
(
    const scalarField& values
) const
{
    const fvMesh& mesh = reconstruction_.mesh();
    const List<symmTensor> secondDerivative(cellSecondDerivatives(values));

    List<FixedList<scalar, 10>> result(mesh.nFaces());
    forAll(result, faceI)
    {
        const label owner = mesh.faceOwner()[faceI];
        const label other = otherCell_[faceI];
        const vector delta = mesh.C()[other] - mesh.C()[owner];
        const scalar distance = mag(delta);
        const vector e = delta/distance;
        const symmTensor J =
            (secondDerivative[other] - secondDerivative[owner])/distance;
        result[faceI] = directionalThirdDerivative(J, e);
    }
    return result;
}


void quadraticGradientCorrection::makeFaceJumpResponse()
{
    const fvMesh& mesh = reconstruction_.mesh();
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());

    // Cubic response of a one-sided face-centre value: the reconstructed value
    // at x_f of a monomial centred at x_f, whose exact value there is zero
    auto& ownerValueCoeffs =
        compactListListCRef(reconstruction_.ownerFaceCentreValueCoeffs());
    auto& neighbourValueCoeffs =
        compactListListCRef(reconstruction_.neighbourFaceCentreValueCoeffs());
    auto valueResponse = [&]
    (
        const label cellI,
        const UList<scalar>& valueCoeffs,
        const point& xf,
        const label term
    )
    {
        const UList<label>& cells = stencils[cellI];
        scalar v = valueCoeffs[cells.size()]*sample(cellI, xf, term);
        forAll(cells, j)
        {
            v += valueCoeffs[j]*sample(cells[j], xf, term);
        }
        return v;
    };

    faceJumpResponse_.setSize(mesh.nFaces());
    forAll(faceJumpResponse_, faceI)
    {
        const label owner = mesh.faceOwner()[faceI];
        const point& xf = mesh.faceCentres()[faceI];
        for (label term = 0; term < 10; ++term)
        {
            faceJumpResponse_[faceI][term] =
              - valueResponse(owner, ownerValueCoeffs[faceI], xf, term);
            if (faceI < mesh.nInternalFaces())
            {
                faceJumpResponse_[faceI][term] += valueResponse
                (
                    mesh.neighbour()[faceI], neighbourValueCoeffs[faceI], xf, term
                );
            }
        }
    }
}


void quadraticGradientCorrection::makeCellFitCoeffs()
{
    const fvMesh& mesh = reconstruction_.mesh();
    const vectorField& C = mesh.C();
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());

    // All monomials through degree three; the ten cubic terms come last
    List<FixedList<label, 3>> basis(20, FixedList<label, 3>(label(0)));
    label n = 1;
    for (direction d = 0; d < 3; ++d)
    {
        basis[n++][d] = 1;
    }
    for (label term = 0; term < 6; ++term, ++n)
    {
        for (direction d = 0; d < 3; ++d)
        {
            basis[n][d] = secondPowers[term][d];
        }
    }
    for (label term = 0; term < 10; ++term, ++n)
    {
        for (direction d = 0; d < 3; ++d)
        {
            basis[n][d] = cubicPowers[term][d];
        }
    }

    cellCubicCoeffs_.setSize(mesh.nCells());
    scalar worstCondition = 0;
    forAll(stencils, cellI)
    {
        // Stencil cells first, the cell itself last
        const UList<label>& cells = stencils[cellI];
        const label nS = cells.size() + 1;
        scalar scale = 0;
        forAll(cells, j)
        {
            scale = max(scale, mag(C[cells[j]] - C[cellI]));
        }

        // Unweighted least-squares cubic fit in coordinates scaled by the
        // stencil radius; every term of a sample is homogeneous in length
        Eigen::MatrixXd A(nS, 20);
        for (label j = 0; j < nS; ++j)
        {
            const label cellJ = j < cells.size() ? cells[j] : cellI;
            for (label col = 0; col < 20; ++col)
            {
                const label p[3] = {basis[col][0], basis[col][1], basis[col][2]};
                A(j, col) = sampleMonomial(cellJ, C[cellI], p)
                   /Foam::pow(scale, p[0] + p[1] + p[2]);
            }
        }
        Eigen::JacobiSVD<Eigen::MatrixXd> svd
        (
            A, Eigen::ComputeThinU | Eigen::ComputeThinV
        );
        const Eigen::VectorXd& sv = svd.singularValues();
        if (nS < 20 || sv(19) <= 1e-10*sv(0))
        {
            FatalErrorInFunction
                << "cellFit cubic stencil is rank deficient at cell " << cellI
                << abort(FatalError);
        }
        worstCondition = max(worstCondition, scalar(sv(0)/sv(19)));
        const Eigen::MatrixXd W =
            svd.matrixV()*sv.cwiseInverse().asDiagonal()
           *svd.matrixU().transpose();

        List<FixedList<scalar, 10>>& coeffs = cellCubicCoeffs_[cellI];
        coeffs.setSize(nS);
        const scalar scale3 = scale*scale*scale;
        for (label j = 0; j < nS; ++j)
        {
            for (label m = 0; m < 10; ++m)
            {
                coeffs[j][m] = W(10 + m, j)/scale3;
            }
        }
    }

    Info<< "Curvature correction for p=2: cellFit cubic least squares on the "
        << "reconstruction stencil, maximum condition number = "
        << worstCondition << nl << endl;
}


void quadraticGradientCorrection::makeFacePatchCoeffs()
{
    const fvMesh& mesh = reconstruction_.mesh();
    const vectorField& C = mesh.C();
    auto& stencils = compactListListCRef(reconstruction_.stencil().cellsStencil());
    auto& coeffs = compactListListCRef(reconstruction_.cellSecondGradCoeffs());
    const labelListList& cellCells = mesh.cellCells();

    // Reconstructed second derivatives of each cubic monomial centred at the
    // cell itself. For a monomial centred at x_f, add its exact second
    // derivative at the cell centre; the quadratic remainder is reproduced.
    List<FixedList<symmTensor, 10>> bias(mesh.nCells());
    forAll(stencils, cellI)
    {
        for (label term = 0; term < 10; ++term)
        {
            const scalar ownSample = sample(cellI, C[cellI], term);
            symmTensor s(symmTensor::zero);
            forAll(stencils[cellI], j)
            {
                s += coeffs[cellI][j]
                   *(sample(stencils[cellI][j], C[cellI], term) - ownSample);
            }
            bias[cellI][term] = s;
        }
    }

    patchCells_.setSize(mesh.nFaces());
    patchCoeffs_.setSize(mesh.nFaces());
    patchJumpCoeffs_.setSize(mesh.nFaces());
    scalar worstCondition = 0;

    forAll(patchCoeffs_, faceI)
    {
        const label owner = mesh.faceOwner()[faceI];
        const label other = otherCell_[faceI];
        labelHashSet set;
        set.insert(owner);
        set.insert(other);
        set.insert(cellCells[owner]);
        set.insert(cellCells[other]);
        patchCells_[faceI] = set.sortedToc();
        const labelList& cells = patchCells_[faceI];
        const label nS = cells.size();

        // Response of the second-derivative deviations from their patch mean
        // to the ten unit cubic modes centred at the face centre
        Eigen::MatrixXd M(6*nS, 10);
        for (label term = 0; term < 10; ++term)
        {
            List<symmTensor> H(nS);
            symmTensor mean(symmTensor::zero);
            forAll(cells, k)
            {
                const vector r = C[cells[k]] - mesh.faceCentres()[faceI];
                symmTensor exact(symmTensor::zero);
                label comp = 0;
                for (direction a = 0; a < 3; ++a)
                {
                    for (direction b = a; b < 3; ++b, ++comp)
                    {
                        label p[3] =
                        {
                            cubicPowers[term][0],
                            cubicPowers[term][1],
                            cubicPowers[term][2]
                        };
                        --p[a];
                        --p[b];
                        if (p[a] >= 0 && p[b] >= 0)
                        {
                            exact[comp] = monomial(r, p);
                        }
                    }
                }
                H[k] = bias[cells[k]][term] + exact;
                mean += H[k]/nS;
            }
            forAll(cells, k)
            {
                for (direction d = 0; d < 6; ++d)
                {
                    M(6*k + d, term) = H[k][d] - mean[d];
                }
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd
        (
            M, Eigen::ComputeThinU | Eigen::ComputeThinV
        );
        const Eigen::VectorXd& sv = svd.singularValues();
        if (sv(9) <= 1e-10*sv(0))
        {
            FatalErrorInFunction
                << "facePatch third-derivative response is rank deficient at"
                << " face " << faceI << abort(FatalError);
        }
        worstCondition = max(worstCondition, scalar(sv(0)/sv(9)));

        // Pseudo-inverse: estimated third derivatives from deviations
        const Eigen::MatrixXd W =
            svd.matrixV()*sv.cwiseInverse().asDiagonal()
           *svd.matrixU().transpose();

        // Gradient increment -faceBeta.T; fold the patch mean into the
        // coefficients so they act directly on the cell second derivatives
        List<FixedList<vector, 6>>& G = patchCoeffs_[faceI];
        G.setSize(nS);
        FixedList<vector, 6> meanG(vector::zero);
        forAll(cells, k)
        {
            for (direction d = 0; d < 6; ++d)
            {
                vector g(vector::zero);
                for (label term = 0; term < 10; ++term)
                {
                    g -= faceBeta_[faceI][term]*W(term, 6*k + d);
                }
                G[k][d] = g;
                meanG[d] += g/nS;
            }
        }
        forAll(cells, k)
        {
            for (direction d = 0; d < 6; ++d)
            {
                G[k][d] -= meanG[d];
            }
        }

        const FixedList<scalar, 10>& jumpResponse = faceJumpResponse_[faceI];
        List<FixedList<scalar, 6>>& J = patchJumpCoeffs_[faceI];
        J.setSize(nS);
        FixedList<scalar, 6> meanJ(scalar(0));
        forAll(cells, k)
        {
            for (direction d = 0; d < 6; ++d)
            {
                scalar c = 0;
                for (label term = 0; term < 10; ++term)
                {
                    c += jumpResponse[term]*W(term, 6*k + d);
                }
                J[k][d] = c;
                meanJ[d] += c/nS;
            }
        }
        forAll(cells, k)
        {
            for (direction d = 0; d < 6; ++d)
            {
                J[k][d] -= meanJ[d];
            }
        }
    }

    Info<< "Curvature correction for p=2: facePatch third-derivative fit, "
        << "maximum response condition number = " << worstCondition
        << nl << endl;
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

quadraticGradientCorrection::quadraticGradientCorrection
(
    const leastSquaresScheme& scheme
)
:
    reconstruction_(scheme),
    nTerms_(scheme.polynomialOrder() == 1 ? 6 : 10),
    variationCoeffs_(),
    responseInverse_(scheme.polynomialOrder() == 1 ? scheme.mesh().nCells() : 0),
    momentQuadrature_(scheme.mesh(), scheme.polynomialOrder() + 1, 0, true),
    faceBeta_(scheme.mesh().nFaces()),
    otherCell_(scheme.polynomialOrder() == 2 ? scheme.mesh().nFaces() : 0, -1),
    faceIntegration_(scheme.curvatureCorrectionFaceIntegration()),
    faceMoment_()
{
    const fvMesh& mesh = scheme.mesh();
    if
    (
        Pstream::parRun() || mesh.nGeometricD() != 3
     || (scheme.polynomialOrder() != 1 && scheme.polynomialOrder() != 2)
    )
    {
        FatalErrorInFunction
            << "curvatureCorrection requires serial 3-D, polynomialOrder 1 or 2"
            << abort(FatalError);
    }
    const word& recovery = scheme.curvatureCorrectionRecovery();
    const bool facePatch = recovery == "facePatch";
    const bool cellFit = recovery == "cellFit";
    if ((facePatch || cellFit) && scheme.polynomialOrder() != 2)
    {
        FatalErrorInFunction
            << "curvatureCorrectionRecovery " << recovery
            << " requires polynomialOrder 2" << abort(FatalError);
    }
    if (faceIntegration_ && !cellFit)
    {
        FatalErrorInFunction
            << "curvatureCorrectionFaceIntegration requires "
            << "curvatureCorrectionRecovery cellFit" << abort(FatalError);
    }
    forAll(mesh.boundary(), patchI)
    {
        const word patchType = mesh.boundaryMesh()[patchI].type();
        if
        (
            mesh.boundary()[patchI].coupled()
         || patchType == "symmetry" || patchType == "symmetryPlane"
        )
        {
            FatalErrorInFunction
                << "curvatureCorrection does not yet support " << patchType
                << " patches" << abort(FatalError);
        }
    }

    if (nTerms_ == 6)
    {
        const vectorField& C = mesh.C();
        auto& stencils = compactListListCRef(scheme.stencil().cellsStencil());
        labelList sizes(mesh.nCells());
        forAll(sizes, cellI)
        {
            sizes[cellI] = stencils[cellI].size();
        }
        variationCoeffs_ = CompactListList<vector>(sizes);

        // A linear fit of derivative differences uses the existing stencil.
        forAll(stencils, cellI)
        {
            const auto& cells = stencils[cellI];
            scalar length = 0;
            forAll(cells, j)
            {
                length = max(length, mag(C[cells[j]] - C[cellI]));
            }
            Eigen::MatrixXd A(cells.size(), 3);
            Eigen::VectorXd w(cells.size());
            forAll(cells, j)
            {
                const vector r = (C[cells[j]] - C[cellI])/length;
                w(j) = 1/max(mag(r), SMALL);
                for (direction d = 0; d < 3; ++d)
                {
                    A(j,d) = w(j)*r[d];
                }
            }
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
            if (qr.rank() != 3)
            {
                FatalErrorInFunction << "Curvature stencil has rank < 3 at cell "
                    << cellI << abort(FatalError);
            }
            const Eigen::MatrixXd inverse = qr.solve
            (
                Eigen::MatrixXd::Identity(cells.size(), cells.size())
            );
            forAll(cells, j)
            {
                variationCoeffs_[cellI][j] =
                    vector(inverse(0,j), inverse(1,j), inverse(2,j))*w(j)/length;
            }
            responseInverse_[cellI] = scalarSquareMatrix(nTerms_, scalar(0));
        }

        // Account for the first-omitted-degree bias of reconstructed derivatives.
        makeResponse(scheme.cellGradCoeffs());
        scalar worstCondition = 0;
        forAll(stencils, cellI)
        {
            Eigen::MatrixXd response(nTerms_,nTerms_);
            for (label i = 0; i < nTerms_; ++i)
            {
                for (label j = 0; j < nTerms_; ++j)
                {
                    response(i,j) = responseInverse_[cellI](i,j);
                }
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd
            (
                response, Eigen::ComputeFullU | Eigen::ComputeFullV
            );
            svd.setThreshold(1e-10);
            if (svd.rank() != nTerms_)
            {
                FatalErrorInFunction << "Derivative response is rank deficient at cell "
                    << cellI << abort(FatalError);
            }
            worstCondition = max
            (
                worstCondition,
                scalar(svd.singularValues()(0)/svd.singularValues()(nTerms_ - 1))
            );
            const Eigen::MatrixXd inverse = svd.solve(Eigen::MatrixXd::Identity(nTerms_,nTerms_));
            for (label i = 0; i < nTerms_; ++i)
            {
                for (label j = 0; j < nTerms_; ++j)
                {
                    responseInverse_[cellI](i,j) = inverse(i,j);
                }
            }
        }

        Info<< "Curvature correction for p=1: maximum derivative-response "
            << "condition number = " << worstCondition << nl << endl;
    }
    else
    {
        // Use the actual adjacent cell internally. At a boundary, choose
        // the face-connected cell most aligned with the inward direction.
        forAll(otherCell_, faceI)
        {
            const label owner = mesh.faceOwner()[faceI];
            if (faceI < mesh.nInternalFaces())
            {
                otherCell_[faceI] = mesh.neighbour()[faceI];
            }
            else
            {
                const vector inward = -mesh.faceAreas()[faceI];
                scalar bestAlignment = 0;
                const labelList& neighbours = mesh.cellCells()[owner];
                forAll(neighbours, j)
                {
                    const label cellI = neighbours[j];
                    const vector delta = mesh.C()[cellI] - mesh.C()[owner];
                    const scalar alignment = (delta & inward)/max(mag(delta), VSMALL);
                    if (alignment > bestAlignment)
                    {
                        bestAlignment = alignment;
                        otherCell_[faceI] = cellI;
                    }
                }
                if (otherCell_[faceI] < 0)
                {
                    // A skewed boundary cell can have no inward immediate
                    // neighbour. Still use only a pair: select the nearest
                    // inward centre from its existing reconstruction stencil.
                    const auto& cells = compactListListCRef
                    (
                        scheme.stencil().cellsStencil()
                    )[owner];
                    scalar nearest = GREAT;
                    forAll(cells, j)
                    {
                        const vector delta = mesh.C()[cells[j]] - mesh.C()[owner];
                        if ((delta & inward) > 0 && mag(delta) < nearest)
                        {
                            nearest = mag(delta);
                            otherCell_[faceI] = cells[j];
                        }
                    }
                }
            }
            if
            (
                otherCell_[faceI] < 0
             || mag(mesh.C()[otherCell_[faceI]] - mesh.C()[owner]) <= VSMALL
            )
            {
                FatalErrorInFunction
                    << "No noncoincident inward cell pair for face " << faceI
                    << abort(FatalError);
            }
        }
        if (!facePatch && !cellFit)
        {
            Info<< "Curvature correction for p=2: direct two-cell directional "
                << "second-derivative differences (no response calibration)"
                << nl << endl;
        }
    }

    auto& faceStencils = compactListListCRef(scheme.faceGradStencil());
    auto& faceCoeffs = const_cast<List<CompactListList<vector>>&>(scheme.faceGradCoeffs());
    auto& points = compactListListCRef(scheme.quadrature().faceQuadPoints());
    const kExactLeastSquares* kExact = isA<kExactLeastSquares>(scheme)
      ? &refCast<const kExactLeastSquares>(scheme) : nullptr;
    auto& weights = compactListListCRef(scheme.quadrature().faceQuadWeights());

    forAll(points, faceI)
    {
        faceBeta_[faceI] = vector::zero;
        const scalar area = mag(mesh.faceAreas()[faceI]);
        forAll(points[faceI], q)
        {
            for (label term = 0; term < nTerms_; ++term)
            {
                vector beta(vector::zero);
                forAll(faceStencils[faceI], j)
                {
                    beta += faceCoeffs[faceI][q][j]
                       *sample(faceStencils[faceI][j], points[faceI][q], term);
                }
                // MLS boundary samples are at the evaluation point, where
                // these monomials vanish. Include all k-exact boundary data.
                if (kExact)
                {
                    const auto& addresses = kExact->faceBoundaryDataStencil()[faceI];
                    auto& dataCoeffs = const_cast<List<CompactListList<vector>>&>
                    (
                        kExact->faceBoundaryDataCoeffs()
                    );
                    forAll(addresses, j)
                    {
                        const auto& a = addresses[j];
                        beta += dataCoeffs[faceI][q][j]*monomial
                        (
                            points[a[0]][a[1]] - points[faceI][q], powers(term)
                        );
                    }
                }
                faceBeta_[faceI][term] += weights[faceI][q]*beta/area;
            }
        }
    }

    if (scheme.curvatureCorrectionBeta() == "faceCentre")
    {
        // One reconstruction of the face stencil at the face centre, as the
        // alpha stabilisation uses, instead of the quadrature average
        if (!isA<movingLeastSquares>(scheme))
        {
            FatalErrorInFunction
                << "curvatureCorrectionBeta faceCentre is implemented for "
                << "movingLeastSquares only" << abort(FatalError);
        }
        const movingLeastSquares& mls = refCast<const movingLeastSquares>(scheme);
        List<vector> coeffs;
        forAll(faceBeta_, faceI)
        {
            mls.faceGradCoeffsAtPoint(faceI, mesh.faceCentres()[faceI], coeffs);
            if (coeffs.empty())
            {
                continue;
            }
            for (label term = 0; term < nTerms_; ++term)
            {
                vector beta(vector::zero);
                forAll(faceStencils[faceI], j)
                {
                    beta += coeffs[j]
                       *sample(faceStencils[faceI][j], mesh.faceCentres()[faceI], term);
                }
                faceBeta_[faceI][term] = beta;
            }
        }
        Info<< "Curvature correction: beta evaluated at the face centres" << nl;
    }

    if (facePatch || cellFit)
    {
        makeFaceJumpResponse();
    }
    if (facePatch)
    {
        makeFacePatchCoeffs();
    }
    else if (cellFit)
    {
        makeCellFitCoeffs();
    }
    if (faceIntegration_)
    {
        makeFaceMoments();
        Info<< "Curvature correction: face-integration error of the "
            << "quadrature rule is corrected" << nl << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

FixedList<scalar, 10> quadraticGradientCorrection::directionalThirdDerivative
(
    const symmTensor& J,
    const vector& e
)
{
    FixedList<scalar, 10> result;
    const vector Je = J & e;
    const scalar eJe = e & Je;

    // Symmetric completion of T_ijk e_k = J_ij. Its purely transverse
    // components vanish. This is a projection, not a fit of a cubic.
    for (label term = 0; term < 10; ++term)
    {
        direction axes[3];
        label n = 0;
        for (direction axis = 0; axis < 3; ++axis)
        {
            for (label power = 0; power < cubicPowers[term][axis]; ++power)
            {
                axes[n++] = axis;
            }
        }
        const direction i = axes[0], j = axes[1], k = axes[2];
        const tensor fullJ(J);
        result[term] =
            e[i]*fullJ(j, k) + e[j]*fullJ(i, k) + e[k]*fullJ(i, j)
          - e[i]*e[j]*Je[k] - e[i]*e[k]*Je[j] - e[j]*e[k]*Je[i]
          + e[i]*e[j]*e[k]*eJe;
    }
    return result;
}


void quadraticGradientCorrection::addTraction
(
    surfaceVectorField& traction,
    const volVectorField& displacement,
    const dimensionedScalar& mu,
    const dimensionedScalar& lambda,
    const scalar scale,
    const scalar alpha,
    const surfaceScalarField* impKfPtr
) const
{
    if (!(scale >= 0 && scale <= 1))
    {
        FatalErrorInFunction << "Correction scale must be between 0 and 1"
            << abort(FatalError);
    }
    if (!(alpha >= 0) || (alpha > 0 && !impKfPtr))
    {
        FatalErrorInFunction << "Alpha compensation needs alpha >= 0 and impKf"
            << abort(FatalError);
    }
    if (scale == 0)
    {
        return;
    }

    const fvMesh& mesh = reconstruction_.mesh();
    if
    (
        &displacement.mesh() != &mesh || &traction.mesh() != &mesh
     || displacement.dimensions() != dimLength
     || traction.dimensions() != dimPressure
     || mu.dimensions() != dimPressure || lambda.dimensions() != dimPressure
    )
    {
        FatalErrorInFunction << "Incompatible fields or elastic moduli"
            << abort(FatalError);
    }

    tensorField gradientIncrement(mesh.nFaces(), tensor::zero);

    // alphaStab adds alpha*impKf*(uN_f - uP_f)/|n.d| on internal faces and
    // alpha*impKf*(u_b - uP_f)/|n.dL| on fixed-value patches. facePatch and
    // cellFit supply all ten third derivatives; twoCell leaves alpha unchanged.
    vectorField alphaIncrement(mesh.nFaces(), vector::zero);
    scalarField alphaWeight(mesh.nFaces(), 0);
    if (alpha > 0 && !faceJumpResponse_.empty())
    {
        const surfaceScalarField& impKf = *impKfPtr;
        const vectorField& C = mesh.C();
        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            const vector n = mesh.faceAreas()[faceI]/mag(mesh.faceAreas()[faceI]);
            const vector d =
                C[mesh.neighbour()[faceI]] - C[mesh.faceOwner()[faceI]];
            alphaWeight[faceI] = alpha*impKf[faceI]/max(mag(n & d), VSMALL);
        }
        forAll(mesh.boundary(), patchI)
        {
            if (!displacement.boundaryField()[patchI].fixesValue())
            {
                continue;
            }
            const label start = mesh.boundaryMesh()[patchI].start();
            forAll(mesh.boundary()[patchI], i)
            {
                const label faceI = start + i;
                const vector n =
                    mesh.faceAreas()[faceI]/mag(mesh.faceAreas()[faceI]);
                const vector dL =
                    mesh.faceCentres()[faceI] - C[mesh.faceOwner()[faceI]];
                alphaWeight[faceI] = alpha*impKf.boundaryField()[patchI][i]
                   /max(mag(n & dL), VSMALL);
            }
        }
    }

    for (direction componentI = 0; componentI < 3; ++componentI)
    {
        const scalarField values
        (
            displacement.primitiveField().component(componentI)
        );
        vector unit(vector::zero);
        unit[componentI] = 1;

        if (!patchCoeffs_.empty())
        {
            const List<symmTensor> H(cellSecondDerivatives(values));
            forAll(gradientIncrement, faceI)
            {
                const labelList& cells = patchCells_[faceI];
                vector increment(vector::zero);
                scalar jump = 0;
                forAll(cells, k)
                {
                    for (direction d = 0; d < 6; ++d)
                    {
                        increment += patchCoeffs_[faceI][k][d]*H[cells[k]][d];
                        jump += patchJumpCoeffs_[faceI][k][d]*H[cells[k]][d];
                    }
                }
                gradientIncrement[faceI] += increment*unit;
                alphaIncrement[faceI][componentI] -= alphaWeight[faceI]*jump;
            }
            continue;
        }

        if (!cellCubicCoeffs_.empty())
        {
            const List<FixedList<scalar, 10>> T(cellThirdDerivatives(values));
            forAll(gradientIncrement, faceI)
            {
                // Owner and neighbour average; owner only on boundary faces
                FixedList<scalar, 10> Tf(T[mesh.faceOwner()[faceI]]);
                if (faceI < mesh.nInternalFaces())
                {
                    const FixedList<scalar, 10>& Tn = T[mesh.neighbour()[faceI]];
                    for (label m = 0; m < 10; ++m)
                    {
                        Tf[m] = 0.5*(Tf[m] + Tn[m]);
                    }
                }
                vector increment(vector::zero);
                scalar jump = 0;
                for (label m = 0; m < 10; ++m)
                {
                    increment -= faceBeta_[faceI][m]*Tf[m];
                    jump += faceJumpResponse_[faceI][m]*Tf[m];
                }
                if (faceIntegration_)
                {
                    // Integration error of the face-average gradient:
                    // M_jk * d3u/dx_i dx_j dx_k
                    const tensor M(faceMoment_[faceI]);
                    for (direction i = 0; i < 3; ++i)
                    {
                        for (direction j = 0; j < 3; ++j)
                        {
                            for (direction k = 0; k < 3; ++k)
                            {
                                increment[i] += M(j, k)*Tf[cubicTerm(i, j, k)];
                            }
                        }
                    }
                }
                gradientIncrement[faceI] += increment*unit;
                alphaIncrement[faceI][componentI] -= alphaWeight[faceI]*jump;
            }
            continue;
        }

        const auto derivatives = nTerms_ == 6
          ? recoveredDerivatives(values) : faceThirdDerivatives(values);
        forAll(gradientIncrement, faceI)
        {
            const label owner = mesh.faceOwner()[faceI];
            vector increment(vector::zero);
            for (label term = 0; term < nTerms_; ++term)
            {
                scalar derivative = derivatives[nTerms_ == 6 ? owner : faceI][term];
                if (nTerms_ == 6 && faceI < mesh.nInternalFaces())
                {
                    derivative = 0.5*
                    (
                        derivative + derivatives[mesh.neighbour()[faceI]][term]
                    );
                }
                increment -= faceBeta_[faceI][term]*derivative;
            }
            gradientIncrement[faceI] += increment*unit;
        }
    }

    forAll(gradientIncrement, faceI)
    {
        const tensor& dGrad = gradientIncrement[faceI];
        const symmTensor dSigma =
            2*mu.value()*symm(dGrad)
          + lambda.value()*tr(dGrad)*symmTensor::I;
        const vector normal =
            mesh.faceAreas()[faceI]/mag(mesh.faceAreas()[faceI]);
        const vector addition =
            scale*((normal & dSigma) + alphaIncrement[faceI]);
        if (faceI < mesh.nInternalFaces())
        {
            primitiveFieldRef(traction)[faceI] += addition;
        }
        else
        {
            const label patchI = mesh.boundaryMesh().whichPatch(faceI);
            boundaryFieldRef(traction)[patchI]
                [faceI - mesh.boundaryMesh()[patchI].start()] += addition;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
