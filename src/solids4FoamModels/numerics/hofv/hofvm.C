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

#include "hofvm.H"
#include "multiplyCoeff.H"
#include "multiplyCoeffExtended.H"
#include "sparseMatrixTools.H"
#include "cellPointLeastSquaresVectors.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::hofvm::laplacian
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    const scalarField& diffusivity,
    const higherOrderGrad& LRE,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): start" << endl;
    }


    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): end" << endl;
    }
}

void Foam::hofvm::laplacianTranspose
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    //const scalarField& diffusivity,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): start" << endl;
    }


    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): end" << endl;
    }
}


void Foam::hofvm::laplacianTrace
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    //const scalarField& diffusivity,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): start" << endl;
    }


    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): end" << endl;
    }
}

// ************************************************************************* //
