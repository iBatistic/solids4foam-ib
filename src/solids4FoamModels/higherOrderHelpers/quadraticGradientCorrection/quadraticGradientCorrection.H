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

Class
    Foam::quadraticGradientCorrection

Description
    Experimental serial 3-D correction of the first omitted polynomial
    degree: quadratic errors for p=1, cubic errors for p=2. The missing
    For p=1, derivatives are recovered with a geometry-only response
    correction. For p=2, twoCell uses the directional difference of the
    owner and neighbour second derivatives, without another fit or
    calibration. facePatch fits all ten third derivatives to the second
    derivatives of the pair and the face neighbours of both cells, using
    their geometry-only response to cubic monomials, so the correction is
    exact for cubic fields. cellFit fits a cubic to each cell's own
    reconstruction stencil and averages owner and neighbour at the face; it
    is also exact for cubic fields but keeps the residual on the p=2
    stencil. With facePatch or cellFit the correction also removes the
    cubic part of the alpha face-value jump.

    Provides one additional traction per face for constant-property linear
    elasticity. Does not change fGrad, face values, alpha stabilisation or
    reconstruction quadrature. Requires polynomialOrder 1 or 2 and a
    matrix-free Jacobian.

SourceFiles
    quadraticGradientCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticGradientCorrection_H
#define quadraticGradientCorrection_H

#include "leastSquaresScheme.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class quadraticGradientCorrection Declaration
\*---------------------------------------------------------------------------*/

class quadraticGradientCorrection
{
    // Private data

        //- Reconstruction whose leading error is corrected
        const leastSquaresScheme& reconstruction_;

        //- Six second derivatives for p=1, ten third derivatives for p=2
        const label nTerms_;

        //- Linear least-squares weights for spatial changes in curvature
        CompactListList<vector> variationCoeffs_;

        //- Inverse response of the derivative-variation operator
        List<scalarSquareMatrix> responseInverse_;

        //- Auxiliary cell moments; does not replace reconstruction quadrature
        fvMeshQuadrature momentQuadrature_;

        //- Area-averaged omitted-degree gradient response, one row per face
        List<FixedList<vector, 10>> faceBeta_;

        //- Second cell of the directional pair (inward cell on boundaries)
        labelList otherCell_;

        //- facePatch: cells whose second derivatives enter each face
        labelListList patchCells_;

        //- facePatch: gradient increment per second-derivative component
        List<List<FixedList<vector, 6>>> patchCoeffs_;

        //- facePatch: cubic part of the alpha face-value jump per
        //  second-derivative component (neighbour minus owner; minus owner
        //  on boundary faces)
        List<List<FixedList<scalar, 6>>> patchJumpCoeffs_;

        //- facePatch and cellFit: cubic response of the alpha face-value
        //  jump per unit third derivative
        List<FixedList<scalar, 10>> faceJumpResponse_;

        //- cellFit: third-derivative coefficients per cell, stencil cells
        //  first and the cell itself last
        List<List<FixedList<scalar, 10>>> cellCubicCoeffs_;

    // Private Member Functions

        //- Evaluate a factorial-scaled monomial
        static scalar monomial(const vector& r, const label* powers);

        //- Exponents of the omitted-degree monomial for this term
        const label* powers(const label term) const;

        //- Sample a monomial, respecting point or average unknowns
        scalar sample(const label cellI, const point& origin, label term) const;

        //- Sample any monomial through degree three
        scalar sampleMonomial
        (
            const label cellI,
            const point& origin,
            const label* powers
        ) const;

        //- Symmetric spatial derivative of a first-derivative field
        FixedList<scalar, 10> variation
        (
            const UList<vector>& firstDerivative,
            const label cellI
        ) const;

        //- Construct the response to the first omitted polynomial degree
        template<class Type>
        void makeResponse(const CompactListList<Type>& coefficients);

        //- Differentiate reconstructed cell derivatives
        template<class Type>
        List<FixedList<scalar, 10>> rawDerivatives
        (
            const CompactListList<Type>& coefficients,
            const scalarField& values
        ) const;

        //- Recover second derivatives for p=1
        List<FixedList<scalar, 10>> recoveredDerivatives
        (
            const scalarField& values
        ) const;

        //- Reconstructed cell second derivatives used by the p=2 recovery
        List<symmTensor> cellSecondDerivatives(const scalarField& values) const;

        //- Directional third derivatives from two reconstructed curvatures
        List<FixedList<scalar, 10>> faceThirdDerivatives
        (
            const scalarField& values
        ) const;

        //- Cubic response of the alpha face-value jump, one row per face
        void makeFaceJumpResponse();

        //- Construct the facePatch coefficients; requires faceBeta_ and
        //  faceJumpResponse_
        void makeFacePatchCoeffs();

        //- Construct the cellFit coefficients: an unweighted cubic
        //  least-squares fit on each cell's reconstruction stencil
        void makeCellFitCoeffs();

        //- Third derivatives of every cell from the cellFit coefficients
        List<FixedList<scalar, 10>> cellThirdDerivatives
        (
            const scalarField& values
        ) const;

        //- Disallow copying
        quadraticGradientCorrection(const quadraticGradientCorrection&) = delete;
        void operator=(const quadraticGradientCorrection&) = delete;

public:

    // Constructors

        //- Construct and cache geometry-only coefficients
        explicit quadraticGradientCorrection(const leastSquaresScheme& scheme);

    // Member Functions

        //- Complete the six directional derivatives as a symmetric tensor.
        //  e must be a unit vector; purely transverse components are zero.
        static FixedList<scalar, 10> directionalThirdDerivative
        (
            const symmTensor& J,
            const vector& e
        );

        //- Add one elastic traction correction per face. Alpha is untouched.
        //  With facePatch and alpha > 0, also remove the cubic part of the
        //  alphaStab high-order face-value jump, alpha*impKf*jump/|n.d|.
        //  Caller must enforce prescribed traction conditions afterwards.
        void addTraction
        (
            surfaceVectorField& traction,
            const volVectorField& displacement,
            const dimensionedScalar& mu,
            const dimensionedScalar& lambda,
            const scalar scale,
            const scalar alpha = 0,
            const surfaceScalarField* impKfPtr = nullptr
        ) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
