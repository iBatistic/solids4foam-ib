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

#include "triQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Map<triQuadrature::quadratureRule> triQuadrature::rules_;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// Initialize and return quadrature rules
void triQuadrature::constructRules()
{
    if (rules_.size() != 0)
    {
        FatalErrorInFunction
            << "attempt to re-construct rules when they already exist"
            << exit(FatalError);
    }

    // 1 point quadrature, exact for polynomials up to 1 order
    rules_.insert
    (
        1,
        quadratureRule
        {
            List<point>{vector(1.0/3.0, 1.0/3.0, 1.0/3.0)},
            List<scalar>{1.0}
        }
    );

    // 3 point quadrature, exact for polynomials up to 2 order
    rules_.insert
    (
        2,
        quadratureRule
        {
            List<point>
            {
                point(2.0/3.0, 1.0/6.0, 1.0/6.0),
                point(1.0/6.0, 2.0/3.0, 1.0/6.0),
                point(1.0/6.0, 1.0/6.0, 2.0/3.0)
            },
            List<scalar>{1.0/3.0, 1.0/3.0, 1.0/3.0}
        }
    );

    // 4 point quadrature, exact for polynomials up to 3 order
    rules_.insert
    (
        3,
        quadratureRule
        {
            List<point>
            {
                point(1.0/3.0, 1.0/3.0, 1.0/3.0),
                point(0.6, 0.2, 0.2),
                point(0.2, 0.6, 0.2),
                point(0.2, 0.2, 0.6)
            },
            List<scalar>{-9.0/16.0, 25.0/48.0, 25.0/48.0, 25.0/48.0}
        }
    );

    // 6 point quadrature, exact for polynomials up to 4 order
    rules_.insert
    (
        4,
        quadratureRule
        {
            List<point>
            {
                point(0.108103018168070, 0.445948490915965, 0.445948490915965),
                point(0.445948490915965, 0.108103018168070, 0.445948490915965),
                point(0.445948490915965, 0.445948490915965, 0.108103018168070),
                point(0.816847572980459, 0.091576213509771, 0.091576213509771),
                point(0.091576213509771, 0.816847572980459, 0.091576213509771),
                point(0.091576213509771, 0.091576213509771, 0.816847572980459)
            },
            List<scalar>
            {
                0.223381589678011,
                0.223381589678011,
                0.223381589678011,
                0.109951743655322,
                0.109951743655322,
                0.109951743655322
            }
        }
    );

    // 7 point quadrature, exact for polynomials up to 5 order
    rules_.insert
    (
        5,
        quadratureRule
        {
            List<point>
            {
                point(1.0/3.0, 1.0/3.0, 1.0/3.0),
                point(0.059715871789770, 0.470142064105115, 0.470142064105115),
                point(0.470142064105115, 0.059715871789770, 0.470142064105115),
                point(0.470142064105115, 0.470142064105115, 0.059715871789770),
                point(0.797426985353087, 0.101286507323456, 0.101286507323456),
                point(0.101286507323456, 0.797426985353087, 0.101286507323456),
                point(0.101286507323456, 0.101286507323456, 0.797426985353087),
            },
            List<scalar>
            {
                 0.225000000000000,
                 0.132394152788506,
                 0.132394152788506,
                 0.132394152788506,
                 0.125939180544827,
                 0.125939180544827,
                 0.125939180544827
            }
        }
    );

    // 12 point quadrature, exact for polynomials up to 6 order
    rules_.insert
    (
        6,
        quadratureRule
        {
            List<point>
            {
                point(0.501426509658179, 0.249286745170910, 0.249286745170910),
                point(0.249286745170910, 0.501426509658179, 0.249286745170910),
                point(0.249286745170910, 0.249286745170910, 0.501426509658179),
                point(0.873821971016996, 0.063089014491502, 0.063089014491502),
                point(0.063089014491502, 0.873821971016996, 0.063089014491502),
                point(0.063089014491502, 0.063089014491502, 0.873821971016996),
                point(0.053145049844817, 0.310352451033784, 0.636502499121399),
                point(0.053145049844817, 0.636502499121399, 0.310352451033784),
                point(0.310352451033784, 0.053145049844817, 0.636502499121399),
                point(0.636502499121399, 0.053145049844817, 0.310352451033784),
                point(0.310352451033784, 0.636502499121399, 0.053145049844817),
                point(0.636502499121399, 0.310352451033784, 0.053145049844817)
            },
            List<scalar>
            {
                0.116786275726379,
                0.116786275726379,
                0.116786275726379,
                0.050844906370207,
                0.050844906370207,
                0.050844906370207,
                0.082851075618374,
                0.082851075618374,
                0.082851075618374,
                0.082851075618374,
                0.082851075618374,
                0.082851075618374
            }
        }
    );
}


const Map<triQuadrature::quadratureRule>& triQuadrature::rules()
{
    if (rules_.size() == 0)
    {
        constructRules();
    }

    return rules_;
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

tmp<Field<point>> triQuadrature::barycentricToPoint
(
    const List<point>& localPts
) const
{
    tmp<Field<point>> tglobalPts(new Field<point>(localPts.size()));
    Field<point>& globalPts = tglobalPts.ref();

    forAll(globalPts, pointI)
    {
        globalPts[pointI] =
            localPts[pointI].x()*this->a()
          + localPts[pointI].y()*this->b()
          + localPts[pointI].z()*this->c();
    }

    return tglobalPts;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


triQuadrature::triQuadrature
(
    const triPoints& pts,
    const label& order
)
:
    quadrature(order),
    triPoints(pts)
{
    // Check if the requested integration order exists in the rules
    if (!rules().found(order))
    {
	FatalErrorInFunction
            << "Quadrature for " << order << " order not implemented"
            << abort(FatalError);
    }

    weights_ = rules()[order].weights;
    points_ = barycentricToPoint(rules()[order].points);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const List<point>& triQuadrature::points() const
{
    return points_;
}


const List<point> triQuadrature::points()
{
    return points_;
}


const List<scalar>& triQuadrature::weights() const
{
    return weights_;
}


const List<scalar> triQuadrature::weights()
{
     return weights_;
}


label triQuadrature::nPoints() const
{
    return points_.size();
}


label triQuadrature::nPoints(label order)
{
    // Check if the requested integration order exists in the rules
    if (!rules().found(order))
    {
	FatalErrorInFunction
            << "Quadrature for " << order << " order not implemented"
            << abort(FatalError);
    }
    return rules()[order].points.size();
}

} // End namespace Foam

// ************************************************************************* //
