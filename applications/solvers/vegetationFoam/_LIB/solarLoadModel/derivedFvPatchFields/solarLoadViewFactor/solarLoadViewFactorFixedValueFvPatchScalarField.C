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

\*---------------------------------------------------------------------------*/

#include "solarLoadViewFactorFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    solarRadiationCoupledBase(patch(), "undefined", scalarField::null()),
    Qso_(p.size(), 0.0)
{}


Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const solarLoadViewFactorFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    solarRadiationCoupledBase
    (
        patch(),
        ptf.albedoMethod(),
        ptf.albedo_
    ),
    Qso_(ptf.Qso_)
{}


Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    solarRadiationCoupledBase(p, dict),
    Qso_("Qso", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );

    }
    else
    {
         fvPatchScalarField::operator=(0.0);
    }
}


Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const solarLoadViewFactorFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    solarRadiationCoupledBase
    (
        ptf.patch(),
        ptf.albedoMethod(),
        ptf.albedo_
    ),
    Qso_(ptf.Qso_)
{}


Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const solarLoadViewFactorFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    solarRadiationCoupledBase
    (
        ptf.patch(),
        ptf.albedoMethod(),
        ptf.albedo_
    ),
    Qso_(ptf.Qso_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Do nothing

    if (debug)
    {
        scalar Q = gSum((*this)*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << " heat transfer rate:" << Q
            << " wall radiative heat flux "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    solarRadiationCoupledBase::write(os);
    Qso_.writeEntry("Qso", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solarLoad
{
    makePatchTypeField
    (
        fvPatchScalarField,
        solarLoadViewFactorFixedValueFvPatchScalarField
    );
}
}


// ************************************************************************* //
