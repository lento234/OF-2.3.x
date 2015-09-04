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

#include "vegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(vegetationModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
vegetationModel::vegetationModel
(
    const volVectorField& U,
    const volScalarField& a
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
    C_
    (
        vegetationProperties_.lookup("C")
    ),
    Cdf_
    (
        vegetationProperties_.lookup("Cdf")
    ),
    l_
    (
        vegetationProperties_.lookup("l")
    ),
    Tl_
    (
        vegetationProperties_.lookup("Tl")
    ),
    UMin_("UMin", dimVelocity, SMALL),
    a_(a)
    {
        Info << "Defined custom vegetation model" << endl;
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// return momentum source term
tmp<fvVectorMatrix> vegetationModel::Su(volVectorField& U) const
{
    return
    (
        - fvm::SuSp(0.5*Cdf_*a_*mag(U), U)
    );
}

tmp<fvScalarMatrix> vegetationModel::Sh(volVectorField& U, volScalarField& T) const
{
    // Calculate magnitude of velocity and bounding above Umin
    volScalarField magU("magU", mag(U));
    bound(magU, UMin_);

    // Calculating aerodynamic boundary layer resistance
    volScalarField ra("ra", C_*pow(l_/magU, 0.5));

    Info << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WIP !!!!!!!!!!!!!!!!!!!!!!!!" << endl;

    return
    (
        //fvm::SuSp((2*a_*(T-Tl_)/ra)/T,T)
        fvm::SuSp((2*a_*(Tl_-T)/ra)/T,T)
    );
}

// return temperature source term

bool vegetationModel::read()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam

// ************************************************************************* //
