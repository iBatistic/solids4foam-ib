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

Description
    Implicit functions for cell-centre high order finite volume calculations:
    High Order Finite Volume Matrix (hofvm).

    Note that these functions only set the matrix coefficients and they do not
    modify the right hand side.

SourceFiles
    hofvm.C

Author
    Ivan Batistic & Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef hofvm_H
#define hofvm_H

#include "volFields.H"
#include "pointFields.H"
#include "sparseMatrix.H"
#include "scalarMatrices.H"
#include "RectangularMatrix.H"
#include "cellPointLeastSquaresVectors.H"
#include "SparseMatrixTemplate.H"
#include "higherOrderGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace hofvm
{
    // Add coefficients to the matrix for the Laplacian of a vector field with
    // a scalar diffusivity field defined at cell centres
    void laplacian
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        const scalarField& diffusivity,
	const higherOrderGrad& LRE,
        const bool debug = false
    );

    // Add coefficients to the matrix for the Laplacian transpose
    void laplacianTranspose
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        //const scalarField& diffusivity,
        const bool debug = false
    );

    // Add coefficients to the matrix for the Laplacian trace
    void laplacianTrace
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        //const scalarField& diffusivity,
        const bool debug = false
    );

} // End namespace hofvm

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
