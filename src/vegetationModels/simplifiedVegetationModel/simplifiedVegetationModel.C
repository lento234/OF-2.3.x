/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "simplifiedVegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(simplifiedVegetationModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
simplifiedVegetationModel::simplifiedVegetationModel
(
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& q
):
    IOdictionary
    (
        IOobject
        (
            "vegetationProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    vegetationProperties_(*this),
    runTime_(U.time()),
    mesh_(U.mesh()),
    cpa_
    (
        vegetationProperties_.lookup("cpa")
    ),
    E_
    (
        vegetationProperties_.lookup("E")
    ),
    Qs_
    (
        vegetationProperties_.lookup("Qs")
    ),
    rhoa_
    (
        vegetationProperties_.lookup("rhoa")
    ),
    Cf_
    (
        IOobject
        (
            "Cf",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Sh_
    (
        IOobject
        (
            "Sh",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,-1,1,0,0,0), 0.0)
    ),
    Sq_
    (
        IOobject
        (
            "Sq",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    ),
    Su_
    (
        IOobject
        (
            "Su",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    )
    {

        Info << " Defined custom vegetation model: simplified vegetation model" << endl;
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// solve vegetation model
void simplifiedVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& q)
{

    Info << "Vegetation model:  Nothing solved. " << endl;

}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> simplifiedVegetationModel::Sh()
{
    Sh_ = Qs_/(rhoa_*cpa_);
    Sh_.correctBoundaryConditions();
    return Sh_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> simplifiedVegetationModel::Su(volVectorField& U)
{
    return fvm::SuSp(-Cf_*mag(U), U);
}

// return specific humidity source term
tmp<volScalarField> simplifiedVegetationModel::Sq()
{
    Sq_ = E_/rhoa_;
    Sq_.correctBoundaryConditions();
    return Sq_;
}

// -----------------------------------------------------------------------------

bool simplifiedVegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam
