/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "rampedParabolicInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
//#include "fvcMeshPhi.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::rampedParabolicInletFvPatchVectorField::currentScale() const
{
    const Foam::scalar t = this->db().time().timeOutputValue();

    if ( t < rampT_ )
    {
        return (1.0 - cos(constant::mathematical::pi/rampT_*t))/2.0;
    }
    else
    {
        return 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rampedParabolicInletFvPatchVectorField::
rampedParabolicInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ubar_(1.0),
    Ufactor_(1.5),
    H_(1.0),
    fluxCorrection_(false),
    rampT_(1.0)
{}


Foam::rampedParabolicInletFvPatchVectorField::
rampedParabolicInletFvPatchVectorField
(
    const rampedParabolicInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Ubar_(ptf.Ubar_),
    Ufactor_(ptf.Ufactor_),
    H_(ptf.H_),
    fluxCorrection_(ptf.fluxCorrection_),
    rampT_(ptf.rampT_)
{}


Foam::rampedParabolicInletFvPatchVectorField::
rampedParabolicInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ubar_(dict.lookupOrDefault<scalar>("Ubar",1.0)),
    Ufactor_(dict.lookupOrDefault<scalar>("Ufactor",1.5)),
    H_(dict.lookupOrDefault<scalar>("H",1.0)),
    fluxCorrection_(dict.lookupOrDefault<bool>("fluxCorrection",false)),
    rampT_(dict.lookupOrDefault<scalar>("ramp",1.0))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::rampedParabolicInletFvPatchVectorField::
rampedParabolicInletFvPatchVectorField
(
    const rampedParabolicInletFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf),
    Ubar_(pvf.Ubar_),
    Ufactor_(pvf.Ufactor_),
    H_(pvf.H_),
    fluxCorrection_(pvf.fluxCorrection_),
    rampT_(pvf.rampT_)
{}


Foam::rampedParabolicInletFvPatchVectorField::
rampedParabolicInletFvPatchVectorField
(
    const rampedParabolicInletFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    Ubar_(pvf.Ubar_),
    Ufactor_(pvf.Ufactor_),
    H_(pvf.H_),
    fluxCorrection_(pvf.fluxCorrection_),
    rampT_(pvf.rampT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rampedParabolicInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();
    const vectorField xc = pp.faceCentres();

    scalar ramp = currentScale();
    //Info<< "Parabolic inlet ramping factor = " << ramp << endl;

    if (fluxCorrection_) // set average velocity over face
    { 
        const labelListList& faceToEdgeIndex = pp.faceEdges();
        const edgeList& edges = pp.edges();
        const pointField& pts = pp.localPoints();
        forAll(*this, i)
        {
            //--faceToEdgeIndex returns a list of edges for each face i
            scalar ymin(9e9);
            scalar ymax(-9e9);
            forAll(faceToEdgeIndex[i],iedge)
            {
                //Info<< "  edge " << iedge << " : " << edges[faceToEdgeIndex[i][iedge]] << endl;
                //--e is the two local points associated with the current edge iedge
                edge e = edges[faceToEdgeIndex[i][iedge]];
                //Info<< "  edge " << iedge << " " << e 
                //    << " : " << pts[e[0]]
                //    << " , " << pts[e[1]] << endl;
                if (pts[e[0]].component(1) != pts[e[1]].component(1))
                {
                    ymin = min( pts[e[0]].component(1), pts[e[1]].component(1) );
                    ymax = max( pts[e[0]].component(1), pts[e[1]].component(1) );
                    break;
                }

            }
            //Info<< "face " << i << " ymin/max : " << ymin << " " << ymax << endl;

            //this->operator[](i).component(0) = ramp * Ufactor_*Ubar_ * y*(H_-y) / (H_*H_/4);
            this->operator[](i).component(0) = 
                ( H_/2*(ymax*ymax - ymin*ymin) - 1.0/3.0*(Foam::pow3(ymax) - Foam::pow3(ymin)) )
                / (ymax-ymin);
            this->operator[](i).component(0) *= ramp * Ufactor_*Ubar_ / (H_*H_/4);
        }
    }
    else // set velocity based on face center
    {
        forAll(*this, i)
        {
            scalar y( xc[i].component(1) );
            //this->operator[](i).component(0) = Ufactor_*Ubar_ * y*(H_-y) / (H_*H_/4);
            this->operator[](i).component(0) = ramp * Ufactor_*Ubar_ * y*(H_-y) / (H_*H_/4);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::rampedParabolicInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Ubar") << Ubar_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ufactor") << Ufactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("H") << H_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluxCorrection") << fluxCorrection_ << token::END_STATEMENT << nl;
    os.writeKeyword("ramp") << rampT_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rampedParabolicInletFvPatchVectorField
    );
}

// ************************************************************************* //
