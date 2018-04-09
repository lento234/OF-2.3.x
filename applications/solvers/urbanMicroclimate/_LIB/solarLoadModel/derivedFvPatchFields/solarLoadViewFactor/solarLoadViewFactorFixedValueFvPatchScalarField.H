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
    Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiative
    heat flux, \c Qr, for use with the view factor model

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Qro          | external radiative heat flux | yes    |
        emissivityMode | emissivity mode: solidThermo or lookup | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            greyDiffusiveRadiationViewFactor;
        Qro             uniform 0;
        emissivityMode  solidThermo;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::radiation::radiationModel
    Foam::radiation::viewFactor
    Foam::radiationCoupledBase
    Foam::fixedValueFvPatchField

SourceFiles
    greyDiffusiveViewFactorFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solarLoadViewFactorFixedValueFvPatchScalarField_H
#define solarLoadViewFactorFixedValueFvPatchScalarField_H

#include "solarRadiationCoupledBase.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solarLoad
{
/*---------------------------------------------------------------------------*\
      Class solarLoadViewFactorFixedValueFvPatchScalarField declaration
\*---------------------------------------------------------------------------*/

class solarLoadViewFactorFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField,
    public solarRadiationCoupledBase
{
    // Private data

        //- External radiative heat flux
        scalarField Qso_;


public:

    //- Runtime type information
    TypeName("solarLoadRadiationViewFactor");


    // Constructors

        //- Construct from patch and internal field
        solarLoadViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solarLoadViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        // solarLoadViewFactorFixedValueFvPatchScalarField onto a new patch
        solarLoadViewFactorFixedValueFvPatchScalarField
        (
            const solarLoadViewFactorFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solarLoadViewFactorFixedValueFvPatchScalarField
        (
            const solarLoadViewFactorFixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new solarLoadViewFactorFixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solarLoadViewFactorFixedValueFvPatchScalarField
        (
            const solarLoadViewFactorFixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new solarLoadViewFactorFixedValueFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return external radiative heat flux
            const scalarList& Qso()
            {
                return Qso_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
