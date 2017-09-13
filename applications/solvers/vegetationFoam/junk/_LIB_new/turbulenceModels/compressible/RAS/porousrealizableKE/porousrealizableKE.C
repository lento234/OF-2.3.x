/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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


---------------------------------
# Implemented vegetation turbulence
# Added compressible terms (aytac) 07.02.2017
# Modified vegetation terms (Lento)

\*---------------------------------------------------------------------------*/

#include "porousrealizableKE.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(porousrealizableKE, 0);
addToRunTimeSelectionTable(RASModel, porousrealizableKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> porousrealizableKE::rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    volScalarField W
    (
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar("small", dimensionSet(0, 0, -3, 0, 0), SMALL)
        )
    );

    tS.clear();

    volScalarField phis
    (
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
    );
    volScalarField As(sqrt(6.0)*cos(phis));
    volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

    return 1.0/(A0_ + As*Us*k_/(epsilon_+epsilonMin_));
}


tmp<volScalarField> porousrealizableKE::rCmu
(
    const volTensorField& gradU
)
{
    volScalarField S2(2*magSqr(dev(symm(gradU))));
    volScalarField magS(sqrt(S2));
    return rCmu(gradU, S2, magS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

porousrealizableKE::porousrealizableKE
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            coeffDict_,
            4.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.9
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            coeffDict_,
            0.9
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            coeffDict_,
            0.9
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.2
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),
    betaP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaP",
            coeffDict_,
            0.0
        )
    ),
    betaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaD",
            coeffDict_,
            0.0
        )
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    LAD_
    (
        IOobject
        (
            "LAD",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    cd_
    (
        IOobject
        (
            "cd",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Cf_
    (
        IOobject
        (
            "Cf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    ),
    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateMut("mut", mesh_)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    mut_ = rCmu(fvc::grad(U_))*rho_*sqr(k_)/epsilon_;
    mut_.correctBoundaryConditions();

    Cf_ = cd_ * LAD_;
    Cf_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    Info << "Defined custom realizableKE model for porous media" << endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> porousrealizableKE::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - (mut_/rho_)*dev(twoSymm(fvc::grad(U_))),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> porousrealizableKE::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> porousrealizableKE::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(T(fvc::grad(U))))
    );
}


bool porousrealizableKE::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        A0_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        C4_.readIfPresent(coeffDict());
        C5_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        betaP_.readIfPresent(coeffDict());
        betaD_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void porousrealizableKE::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rCmu(fvc::grad(U_))*rho_*sqr(k_)/epsilon_;
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    volScalarField divU(fvc::div(phi_/fvc::interpolate(rho_)));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    volTensorField gradU(fvc::grad(U_));
    volScalarField S2(2*magSqr(dev(symm(gradU))));
    volScalarField magS(sqrt(S2));

    volScalarField eta(magS*k_/epsilon_);
    volScalarField C1(max(eta/(scalar(5) + eta), scalar(0.43)));

    volScalarField G(GName(), mut_*(gradU && dev(twoSymm(gradU))));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*rho_*magS*epsilon_
      - fvm::Sp
        (
            C2_*rho_*epsilon_/(k_ + sqrt((mu()/rho_)*epsilon_)),
            epsilon_
        )
      + fvm::Sp(rho_*betaP_*C4_*Cf_*pow(mag(U_),3)/k_,epsilon_) // TODO: look into this: is it not more stable if positive term is without fvm::Sp (so that it is not in the diagonal)? - ayk
      - fvm::Sp(rho_*betaD_*C5_*Cf_*mag(U_),epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);

    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::SuSp(2.0/3.0*rho_*divU, k_)
      - fvm::Sp(rho_*epsilon_/k_, k_)
      + fvm::Sp(rho_*betaP_*Cf_*pow(mag(U_),3)/k_, k_) // TODO: look into this: is it not more stable if positive term is without fvm::Sp (so that it is not in the diagonal)? - ayk
      - fvm::Sp(rho_*betaD_*Cf_*mag(U_), k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    // Re-calculate viscosity
    mut_ = rCmu(gradU, S2, magS)*rho_*sqr(k_)/(epsilon_+epsilonMin_);
    Info << ">>>>>>>>>>>>>>>>>>>>>>>>>>> survived, dont forget about thiiiis !!!" << endl;
    mut_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
