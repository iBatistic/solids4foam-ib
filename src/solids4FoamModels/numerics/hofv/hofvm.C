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
#include "fvc.H"
#include "multiplyCoeff.H"
#include "multiplyCoeffExtended.H"
#include "sparseMatrixTools.H"
#include "cellPointLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"
#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::hofvm::laplacian
(
    sparseMatrix& matrix,
    vectorField& source,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const higherOrderGrad& LRE,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::hofvm::laplacian(...): start" << endl;
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Diffusion coefficient linearly interpolated to face centres.
    // TO_DO: interpolate diffusivity to quad points using hofvc::interpolate
    const surfaceScalarField gamma(fvc::interpolate(diffusivity));
    const scalarField& gammaI = gamma.internalField();

    // Face quadrature points weights
    const List<List<scalar>>& facesQuadWeights = LRE.faceGaussPointsWeight();

    // Faces stencil
    const List<labelList>& stencils = LRE.globalFaceStencils();

    // Gradient interpolation coefficients
    const List<List<DynamicList<vector>>>& gradCoeffs =
	LRE.QRGradFaceGPCoeffs();

    // Loop over internal faces
    forAll(owner, faceI)
    {
	// Preliminaries
	const vector& faceNormal = n[faceI];

	const scalar gammaMagSf = magSf[faceI] * gammaI[faceI];

	// Face interpolation molecule
	const labelList& faceStencil = stencils[faceI];

	// Face quadrature points weights
	const List<scalar>& faceQuadWeight = facesQuadWeights[faceI];

	// Loop over face quadrature points
	forAll(faceQuadWeight, pointI)
	{
	    // Quad point weight
	    const scalar& quadPointW = faceQuadWeight[pointI];

	    // Loop over interpolation stencil
	    for(label cI = 0; cI < faceStencil.size(); cI++)
	    {
		const label globalCellID = faceStencil[cI];
		const vector& cellGradCoeff = gradCoeffs[faceI][pointI][cI];

		const tensor coeff =
		    I*gammaMagSf*quadPointW*(cellGradCoeff&faceNormal);

		matrix(owner[faceI], globalCellID) += coeff;
		matrix(neighbour[faceI], globalCellID) -= coeff;
	    }
	}
    }

    // Loop over boundary faces
    forAll(mesh.boundaryMesh(), patchI)
    {
	const word& patchType = mesh.boundaryMesh()[patchI].type();
	const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
	const scalarField& pGamma = gamma.boundaryField()[patchI];
	const vectorField& pNormal = mesh.boundary()[patchI].nf();

	const label start = mesh.boundaryMesh()[patchI].start();

	if (patchType == emptyPolyPatch::typeName)
	{
	    // Skip empty patches
	}
	else if (patchType == processorPolyPatch::typeName)
	{
	    NotImplemented;
	}
	else if (patchType == symmetryPolyPatch::typeName)
	{
	    NotImplemented;
	}
	else if
	(
	    isA<solidTractionFvPatchVectorField>(D.boundaryField()[patchI])
	)
        {
	    const solidTractionFvPatchVectorField& tracPatch =
	       refCast<const solidTractionFvPatchVectorField>
               (
	           D.boundaryField()[patchI]
               );

	    const vectorField traction
	    (
	        tracPatch.traction() - pNormal*tracPatch.pressure()
	    );

            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
		// Get global face index
		const label faceID = faceI + start;

		// Face force
		const vector force = pMagSf[faceI] * traction[faceI];

		source[owner[faceID]] -= force;
	    }
        }
	else if
	(
	    isA<fixedDisplacementFvPatchVectorField>(D.boundaryField()[patchI])
	)
	{
            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
		// Preliminaries
		const vector& faceNormal = pNormal[faceI];
		const scalar gammaMagSf = pMagSf[faceI] * pGamma[faceI];

		// Get global face index, needed for lists from LRE class
		const label faceID = faceI + start;

		// Face interpolation molecule
		const labelList& faceStencil = stencils[faceID];

		// Face quadrature points weights
		const List<scalar>& faceQuadWeight = facesQuadWeights[faceID];

		forAll(faceQuadWeight, pointI)
		{
		    // Quad point weight
		    const scalar& quadPointW = faceQuadWeight[pointI];

		    // Loop over interpolation stencil. Last item in faceStencil
		    // is boundary face itself. Treated separately after loop.
		    for(label cI = 0; cI < (faceStencil.size() - 1); cI++)
		    {
			const label globalCellID = faceStencil[cI];
			const vector& cellGradCoeff =
			    gradCoeffs[faceID][pointI][cI];

			const tensor coeff =
			    I*gammaMagSf*quadPointW*(cellGradCoeff&faceNormal);

			matrix(owner[faceID], globalCellID) += coeff;
		    }

		    // Include boundary value to the source
		    {
			const label size = faceStencil.size();

			const vector& cellGradCoeff =
			    gradCoeffs[faceID][pointI][size];

			const tensor coeff =
			    I*gammaMagSf*quadPointW*(cellGradCoeff&faceNormal);

			source[owner[faceID]] -=
			    coeff & D.boundaryField()[patchI][faceI];
		    }
		}
	    }
	}
	else
	{
	    NotImplemented;
	}
    }

    //matrix.print();

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
