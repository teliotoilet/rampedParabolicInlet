/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::rampedParabolicInletFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Implementation started from movingWallVelocityFvPatchVectorField.

    Provides a parabolic inlet profile (according to the parameters specified
    by Turek and Hron 2006. Ramping over time is applied with a 
    (1-cos(pi/2*t))/2 function. 

SeeAlso
    Foam::fixedValueFvPatchVectorField

SourceFiles
    rampedParabolicInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rampedParabolicInletFvPatchVectorField_H
#define rampedParabolicInletFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class rampedParabolicInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

       //- Inlet profile parameters
       scalar Ubar_;        // average velocity
       scalar Ufactor_;     // scales the maximum velocity
       scalar H_;           // channel height

       //- Ramping period
       scalar rampT_;       // inlet profile constant after this period


    // Private Member Functions

        //- Return current scale
        scalar currentScale() const;

public:

    //- Runtime type information
    TypeName("rampedParabolicInlet");


    // Constructors

        //- Construct from patch and internal field
        rampedParabolicInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rampedParabolicInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rampedParabolicInletFvPatchVectorField
        //  onto a new patch
        rampedParabolicInletFvPatchVectorField
        (
            const rampedParabolicInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rampedParabolicInletFvPatchVectorField
        (
            const rampedParabolicInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new rampedParabolicInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rampedParabolicInletFvPatchVectorField
        (
            const rampedParabolicInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new rampedParabolicInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
