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

#include "solarRadiationCoupledBase.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fvPatchFieldMapper.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::solarRadiationCoupledBase::albedoMethodType,
        2
    >::names[] =
    {
        "solidRadiation",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::solarRadiationCoupledBase::albedoMethodType, 2>
    Foam::solarRadiationCoupledBase::albedoMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarRadiationCoupledBase::solarRadiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const scalarField& albedo
)
:
    patch_(patch),
    method_(albedoMethodTypeNames_[calculationType]),
    albedo_(albedo)
{}


Foam::solarRadiationCoupledBase::solarRadiationCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(albedoMethodTypeNames_.read(dict.lookup("albedoMode")))
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            if (!isA<mappedPatchBase>(patch_.patch()))
            {
                FatalIOErrorIn
                (
                    "solarRadiationCoupledBase::solarRadiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    patch type '" << patch_.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << patch_.name()
                    << exit(FatalIOError);
            }

            albedo_ = scalarField(patch_.size(), 0.0);
        }
        break;

        case LOOKUP:
        {
            if (!dict.found("albedo"))
            {
                FatalIOErrorIn
                (
                    "solarRadiationCoupledBase::solarRadiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    albedo key does not exist for patch "
                    << patch_.name()
                    << exit(FatalIOError);
            }
            else
            {
                albedo_ = scalarField("albedo", dict, patch_.size());
            }
        }
        break;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::solarRadiationCoupledBase::albedo() const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_.patch());

            const polyMesh& nbrMesh = mpp.sampleMesh();

            const radiation::radiationModel& radiation =
                nbrMesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );


            const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

            const fvPatch& nbrPatch =
                nbrFvMesh.boundary()[mpp.samplePolyPatch().index()];


            scalarField albedo
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    nbrPatch.index()
                ]
            );
            mpp.distribute(albedo);

            return albedo;

        }
        break;

        case LOOKUP:
        {
            // return local value
            return albedo_;
        }

        default:
        {
            FatalErrorIn
            (
                "solarRadiationCoupledBase::albedo(const scalarField&)"
            )   << "Unimplemented method " << method_ << endl
                << "Please set 'albedo' to one of "
                << albedoMethodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


void Foam::solarRadiationCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("albedoMode") << albedoMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    albedo_.writeEntry("albedo", os);
}


// ************************************************************************* //
