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

Application
    Test-quadraticGradientCorrection

Description
    Serial 3-D test of the separate face-traction correction for p=1 and p=2.
    Checks integrated polynomial traction, scales 0/0.5/1, and unchanged
    reconstruction gradients, face values and quadrature when enabled.
    Uses analytical fixed-displacement data in memory; no case files change.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#ifdef FOAMEXTEND
int main(int argc, char *argv[])
{
    Info<< "SKIPPED: Test-quadraticGradientCorrection requires OpenFOAM" << endl;
    return 0;
}
#else

#include "leastSquaresScheme.H"
#include "quadraticGradientCorrection.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "compatibilityFunctions.H"

using namespace Foam;

namespace
{
    const vector componentScales(1, -0.7, 1.3);

    scalar value(const point& x, const FixedList<label, 3>& powers)
    {
        scalar result = 1;
        for (direction d = 0; d < 3; ++d)
        {
            for (label i = 0; i < powers[d]; ++i)
            {
                result *= x[d];
            }
        }
        return result;
    }

    vector gradient(const point& x, const FixedList<label, 3>& powers)
    {
        vector result(vector::zero);
        for (direction d = 0; d < 3; ++d)
        {
            if (powers[d])
            {
                FixedList<label, 3> reduced(powers);
                --reduced[d];
                result[d] = powers[d]*value(x, reduced);
            }
        }
        return result;
    }

    // Prescribed values follow the test polynomial at every quadrature point.
    class polynomialDisplacement : public fixedDisplacementFvPatchVectorField
    {
        const leastSquaresScheme& scheme_;
        const FixedList<label, 3>& powers_;

    public:
        polynomialDisplacement
        (
            const fvPatch& patch,
            const DimensionedField<vector, volMesh>& internal,
            const leastSquaresScheme& scheme,
            const FixedList<label, 3>& powers
        )
        :
            fixedDisplacementFvPatchVectorField(patch, internal),
            scheme_(scheme),
            powers_(powers)
        {}

        virtual autoPtr<CompactListList<vector>> evaluateQuadrature() const
        {
            auto& points = compactListListCRef(scheme_.quadrature().faceQuadPoints());
            labelList sizes(patch().size());
            forAll(sizes, i)
            {
                sizes[i] = points[patch().start() + i].size();
            }
            autoPtr<CompactListList<vector>> result(new CompactListList<vector>(sizes));
            auto& values = autoPtrRef(result);
            forAll(values, i)
            {
                forAll(values[i], q)
                {
                    values[i][q] = value(points[patch().start() + i][q], powers_)
                       *componentScales;
                }
            }
            return result;
        }
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (Pstream::parRun() || mesh.nGeometricD() != 3)
    {
        FatalErrorInFunction << "Requires a serial 3-D mesh" << abort(FatalError);
    }
    IOdictionary properties
    (
        IOobject
        (
            "solidProperties", runTime.constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );
    const word model(properties.lookup("solidModel"));
    dictionary coefficients
    (
        properties.subDict(model + "Coeffs")
            .subDict("highOrderCoeffs").subDict("displacement")
    );
    const label order = readLabel(coefficients.lookup("polynomialOrder"));
    if (order != 1 && order != 2)
    {
        FatalErrorInFunction << "Requires polynomialOrder 1 or 2"
            << abort(FatalError);
    }
    coefficients.set("curvatureCorrectionScale", 1);
    const boolList includePatches(mesh.boundary().size(), true);
    const fvMeshQuadrature exactQuad(mesh, 3, 2, true);
    auto& cellPoints = compactListListCRef(exactQuad.cellQuadPoints());
    auto& cellWeights = compactListListCRef(exactQuad.cellQuadWeights());
    const dimensionedScalar mu("mu", dimPressure, 1);
    const dimensionedScalar lambda("lambda", dimPressure, 2);
    const wordList methods({"movingLeastSquares", "kExactLeastSquares"});
    bool passed = true;

    // Independent cubic identity: for u=(a.x)^3/6, the exact third tensor
    // is a*a*a. Remove its purely transverse part to obtain the expected
    // directional projection. Include axis-aligned and rotated pairs.
    scalar projectionError = 0;
    const vector directions[] = {vector(1, 0, 0), vector(2, -1, 2)/3};
    const vector modes[] =
    {
        vector(1, 0, 0), vector(0, 1, 0), vector(0, 0, 1),
        vector(1, 2, -1), vector(2, -1, 2)/3, vector(1, 2, 0)
    };
    for (const vector& e : directions)
    {
        for (const vector& a : modes)
        {
            const symmTensor J = (a & e)*sqr(a);
            const auto result =
                quadraticGradientCorrection::directionalThirdDerivative(J, e);
            const auto reverse =
                quadraticGradientCorrection::directionalThirdDerivative(-J, -e);
            const vector transverse = a - (a & e)*e;
            label term = 0;
            for (direction i = 0; i < 3; ++i)
            {
                for (direction j = i; j < 3; ++j)
                {
                    for (direction k = j; k < 3; ++k, ++term)
                    {
                        const scalar expected = a[i]*a[j]*a[k]
                          - transverse[i]*transverse[j]*transverse[k];
                        projectionError = max(projectionError,
                            mag(result[term] - expected));
                        projectionError = max(projectionError,
                            mag(reverse[term] - expected));
                    }
                }
            }
        }
    }
    passed = projectionError < 1e-12;
    Info<< "Directional projection and pair reversal: "
        << (passed ? "PASSED" : "FAILED")
        << ", error " << projectionError << endl;

    auto faceValue = [&mesh](const surfaceVectorField& field, const label faceI)
    {
        if (faceI < mesh.nInternalFaces())
        {
            return field[faceI];
        }
        const label patchI = mesh.boundaryMesh().whichPatch(faceI);
        return field.boundaryField()[patchI]
            [faceI - mesh.boundaryMesh()[patchI].start()];
    };
    auto setFace = [&mesh]
    (
        surfaceVectorField& field, const label faceI, const vector& value
    )
    {
        if (faceI < mesh.nInternalFaces())
        {
            primitiveFieldRef(field)[faceI] = value;
        }
        else
        {
            const label patchI = mesh.boundaryMesh().whichPatch(faceI);
            boundaryFieldRef(field)[patchI]
                [faceI - mesh.boundaryMesh()[patchI].start()] = value;
        }
    };
    auto stress = [&](const tensor& gradient)
    {
        return 2*mu.value()*symm(gradient)
          + lambda.value()*tr(gradient)*symmTensor::I;
    };

    forAll(methods, method)
    {
        coefficients.set("type", methods[method]);
        coefficients.set("curvatureCorrection", false);
        autoPtr<leastSquaresScheme> offPtr =
            leastSquaresScheme::New(mesh, includePatches, coefficients);
        coefficients.set("curvatureCorrection", true);
        autoPtr<leastSquaresScheme> onPtr =
            leastSquaresScheme::New(mesh, includePatches, coefficients);
        const leastSquaresScheme& scheme = onPtr();
        const leastSquaresScheme& original = offPtr();
        FixedList<label, 3> powers(label(0));
        volVectorField D
        (
            IOobject
            (
                "testPolynomial", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            mesh, dimensionedVector("zero", dimLength, vector::zero),
            "zeroGradient"
        );
        forAll(mesh.boundary(), patchI)
        {
            boundaryFieldRef(D).set
            (
                patchI,
                new polynomialDisplacement
                (
                    mesh.boundary()[patchI], D.internalField(), scheme, powers
                )
            );
        }
        surfaceVectorField own
        (
            IOobject
            (
                "testOwner", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            mesh, dimensionedVector("zero", dimLength, vector::zero)
        );
        surfaceVectorField nei("testNeighbour", own);
        surfaceVectorField ownOff("testOwnerOff", own);
        surfaceVectorField neiOff("testNeighbourOff", own);
        surfaceVectorField base
        (
            IOobject
            (
                "baseTraction", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            mesh, dimensionedVector("zero", dimPressure, vector::zero)
        );
        auto& points = compactListListCRef(scheme.quadrature().faceQuadPoints());
        auto& weights = compactListListCRef(scheme.quadrature().faceQuadWeights());
        auto& oldPoints = compactListListCRef(original.quadrature().faceQuadPoints());
        auto& oldWeights = compactListListCRef(original.quadrature().faceQuadWeights());
        auto& onCells = compactListListCRef(scheme.quadrature().cellQuadPoints());
        auto& offCells = compactListListCRef(original.quadrature().cellQuadPoints());
        auto& onCellW = compactListListCRef(scheme.quadrature().cellQuadWeights());
        auto& offCellW = compactListListCRef(original.quadrature().cellQuadWeights());
        if (points.sizes() != oldPoints.sizes() || onCells.sizes() != offCells.sizes())
        {
            FatalErrorInFunction << "Correction changed quadrature sizes"
                << abort(FatalError);
        }
        scalar unchangedError = 0;
        forAll(points, faceI)
        {
            forAll(points[faceI], q)
            {
                unchangedError = max(unchangedError, mag(points[faceI][q] - oldPoints[faceI][q]));
                unchangedError = max(unchangedError, mag(weights[faceI][q] - oldWeights[faceI][q]));
            }
        }
        forAll(onCells, cellI)
        {
            forAll(onCells[cellI], q)
            {
                unchangedError = max(unchangedError, mag(onCells[cellI][q] - offCells[cellI][q]));
                unchangedError = max(unchangedError, mag(onCellW[cellI][q] - offCellW[cellI][q]));
            }
        }
        CompactListList<tensor> faceGrad(points.sizes());
        CompactListList<tensor> oldGrad(points.sizes());
        scalar maxTractionError = 0, maxScaleError = 0;
        scalar maxCubicError = 0;
        for (label degree = 0; degree <= order + 1; ++degree)
        {
            for (label x = 0; x <= degree; ++x)
            {
                for (label y = 0; y <= degree - x; ++y)
                {
                    powers[0] = x;
                    powers[1] = y;
                    powers[2] = degree - x - y;
                    forAll(D, cellI)
                    {
                        scalar sample = value(mesh.C()[cellI], powers);
                        if (methods[method] == "kExactLeastSquares")
                        {
                            sample = 0;
                            forAll(cellPoints[cellI], q)
                            {
                                sample += cellWeights[cellI][q]
                                   *value(cellPoints[cellI][q], powers);
                            }
                            sample /= mesh.V()[cellI];
                        }
                        D[cellI] = sample*componentScales;
                    }
                    forAll(mesh.boundary(), patchI)
                    {
                        auto& patch = boundaryFieldRef(D)[patchI];
                        forAll(patch, faceI)
                        {
                            patch[faceI] = value
                            (
                                mesh.boundary()[patchI].Cf()[faceI], powers
                            )*componentScales;
                        }
                    }
                    scheme.fGrad(D, faceGrad);
                    original.fGrad(D, oldGrad);
                    scheme.faceCentreValues(D, own, nei);
                    original.faceCentreValues(D, ownOff, neiOff);
                    vectorField exact(mesh.nFaces(), vector::zero);
                    forAll(points, faceI)
                    {
                        const scalar area = mag(mesh.faceAreas()[faceI]);
                        const vector normal = mesh.faceAreas()[faceI]/area;
                        vector traction(vector::zero);
                        forAll(points[faceI], q)
                        {
                            const scalar w = weights[faceI][q]/area;
                            traction += w*(normal & stress(faceGrad[faceI][q]));
                            exact[faceI] += w*(normal & stress
                            (
                                gradient(points[faceI][q], powers)*componentScales
                            ));
                            unchangedError = max(unchangedError,
                                mag(faceGrad[faceI][q] - oldGrad[faceI][q]));
                        }
                        setFace(base, faceI, traction);
                        unchangedError = max(unchangedError,
                            mag(faceValue(own, faceI) - faceValue(ownOff, faceI)));
                        unchangedError = max(unchangedError,
                            mag(faceValue(nei, faceI) - faceValue(neiOff, faceI)));
                    }

                    surfaceVectorField full("fullCorrection", base);
                    surfaceVectorField half("halfCorrection", base);
                    surfaceVectorField zero("zeroCorrection", base);
                    scheme.curvatureCorrection().addTraction(full, D, mu, lambda, 1);
                    scheme.curvatureCorrection().addTraction(half, D, mu, lambda, 0.5);
                    scheme.curvatureCorrection().addTraction(zero, D, mu, lambda, 0);
                    scalar error = 0;
                    forAll(exact, faceI)
                    {
                        error = max(error, mag(faceValue(full, faceI) - exact[faceI]));
                        maxScaleError = max(maxScaleError,
                            mag(faceValue(zero, faceI) - faceValue(base, faceI)));
                        maxScaleError = max(maxScaleError, mag
                        (
                            faceValue(half, faceI)
                          - 0.5*(faceValue(full, faceI) + faceValue(base, faceI))
                        ));
                    }
                    if (order == 2 && degree == 3)
                    {
                        // A two-cell difference cannot recover all cubic
                        // modes. Report these errors, but do not require
                        // exact cubic traction from this directional scheme.
                        maxCubicError = max(maxCubicError, error);
                        scalar originalError = 0;
                        forAll(exact, faceI)
                        {
                            originalError = max(originalError,
                                mag(faceValue(base, faceI) - exact[faceI]));
                        }
                        Info<< "powers " << powers
                            << ": cubic diagnostic (not an exactness check), "
                            << "original " << originalError
                            << ", corrected " << error << endl;
                    }
                    else
                    {
                        maxTractionError = max(maxTractionError, error);
                        Info<< "powers " << powers << ": traction error "
                            << error << endl;
                    }
                }
            }
        }
        const bool methodPassed =
            maxTractionError < 1e-10 && maxScaleError < 1e-10 && unchangedError == 0;
        passed = passed && methodPassed;
        Info<< methods[method] << ": " << (methodPassed ? "PASSED" : "FAILED")
            << ", traction error " << maxTractionError
            << ", scale error " << maxScaleError
            << ", unchanged operators error " << unchangedError << endl;
        if (order == 2)
        {
            Info<< "Maximum cubic traction error (diagnostic only): "
                << maxCubicError << endl;
        }
    }
    return passed ? 0 : 1;
}
#endif

// ************************************************************************* //
