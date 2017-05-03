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

#include "modifiedVegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(modifiedVegetationModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
modifiedVegetationModel::modifiedVegetationModel
(
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& w
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
    a1_
    (
        vegetationProperties_.lookup("a1")
    ),
    a2_
    (
        vegetationProperties_.lookup("a2")
    ),
    a3_
    (
        vegetationProperties_.lookup("a3")
    ),
    cpa_
    (
        vegetationProperties_.lookup("cpa")
    ),
    C_
    (
        vegetationProperties_.lookup("C")
    ),
    D0_
    (
        vegetationProperties_.lookup("D0")
    ),
    nEvapSides_
    (
        vegetationProperties_.lookup("nEvapSides")
    ),
    H_
    (
        vegetationProperties_.lookup("H")
    ),
    kc_
    (
        vegetationProperties_.lookup("kc")
    ),
    l_
    (
        vegetationProperties_.lookup("l")
    ),
    Rg0_
    (
        vegetationProperties_.lookup("Rg0")
    ),
    Rl0_
    (
        vegetationProperties_.lookup("Rl0")
    ),
    rhoa_
    (
        vegetationProperties_.lookup("rhoa")
    ),
    rsMin_
    (
        vegetationProperties_.lookup("rsMin")
    ),
    TlMin_("TlMin", dimTemperature, SMALL),
    UMin_("UMin", dimVelocity, SMALL),
    lambda_
    (
        vegetationProperties_.lookup("lambda")
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
    E_
    (
        IOobject
        (
            "E",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),
    ev_
    (
        IOobject
        (
            "ev",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),
    evsat_
    (
        IOobject
        (
            "evsat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
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
    LAI_
    (
        IOobject
        (
            "LAI",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    wsat_
    (
        IOobject
        (
            "wsat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    Ql_
    (
        IOobject
        (
            "Ql",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    Qs_
    (
        IOobject
        (
            "Qs",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    ra_
    (
        IOobject
        (
            "ra",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    rs_
    (
        IOobject
        (
            "rs",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    rhosat_
    (
        IOobject
        (
            "rhosat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    ),
    Rg_
    (
        IOobject
        (
            "Rg",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimensionSet(1,0,-3,0,0,0,0), vector::zero)
    ),
    Rn_
    (
        IOobject
        (
            "Rn",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    ST_
    (
        IOobject
        (
            "ST",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,-1,1,0,0,0), 0.0)
    ),
    Sw_
    (
        IOobject
        (
            "Sw",
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
    ),
    Tl_
    (
        IOobject
        (
            "Tl",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,0,1,0,0,0), 0.0)
    ),
    VPD_
    (
        IOobject
        (
            "VPD",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    )
    {
        // Bounding parameters
        bound(Tl_, TlMin_);


        Info << " Defined custom vegetation model" << endl;

        Cf_ = cd_ * LAD_;

    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double modifiedVegetationModel::calc_evsat(double& T)
{
    // saturated vapor pressure pws - ASHRAE 1.2
    return exp( - 5.8002206e3/T
                + 1.3914993
                - 4.8640239e-2*T
                + 4.1764768e-5*pow(T,2)
                - 1.4452093e-8*pow(T,3)
                + 6.5459673*log(T) );
}

// calc saturated density of water vapour
double modifiedVegetationModel::calc_rhosat(double& T)
{
    return calc_evsat(T)/(461.5*T);
}

// solve radiation
void modifiedVegetationModel::radiation()
{
    // empirical global radiation = f(height) inside and below vegetation
    forAll(Rg_, cellI)
        Rg_[cellI] = vector(0,0, Rg0_.value()*exp(-kc_.value()*LAI_[cellI]));
    Rg_.correctBoundaryConditions();

    volTensorField gradRg = fvc::grad(Rg_);
    gradRg.correctBoundaryConditions();

    // radiation density inside vegetation
    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Rn_[cellI] = gradRg[cellI].component(8) + 0.04*(Rl0_.value()/H_.value()); //gradRg[cellI] && tensor(0,0,0,0,0,0,0,0,1);
    Rn_.correctBoundaryConditions();

}

// solve aerodynamic resistance
void modifiedVegetationModel::resistance(volScalarField& magU, volScalarField& T, volScalarField& w, volScalarField& Tl)
{
    const double p_ = 101325;

    // Calculate magnitude of velocity and bounding above Umin
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            //Aerodynamic resistance
            ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

            // Calculate vapor pressure of air
            ev_[cellI] = p_*w[cellI]/(0.621945+w[cellI]);

            // Calculate sat. vapor pressure of air
            evsat_[cellI] = calc_evsat(T[cellI]);

            // Vapor pressure deficit - kPa
            VPD_[cellI] = evsat_[cellI] - ev_[cellI];

            // Stomatal resistance - type 2
            if ((VPD_[cellI]/1000.0) < D0_.value())
                rs_[cellI] = rsMin_.value();//*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])));
            else
                rs_[cellI] = rsMin_.value();//*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));

        }
    }
    ev_.correctBoundaryConditions();
    evsat_.correctBoundaryConditions();
    VPD_.correctBoundaryConditions();
    ra_.correctBoundaryConditions();
    rs_.correctBoundaryConditions();
}

// solve vegetation model
void modifiedVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& w)
{
    // solve radiation within vegetation
    radiation();

    const double p_ = 101325;

    // Magnitude of velocity
    volScalarField magU("magU", mag(U));

    // Bounding velocity
    bound(magU, UMin_);

    // solve aerodynamic, stomatal resistance
    volScalarField new_Tl("new_Tl", Tl_);

    // info
    Info << "    max leaf temp tl=" << max(T.internalField())
         << "k, iteration i=0" << endl;

    scalar maxError, maxRelError;
    int i;

    // solve leaf temperature, iteratively.
    int maxIter = 500;
    for (i=1; i<=maxIter; i++)
    {
        // Solve aerodynamc, stomatal resistance
        resistance(magU, T, w, new_Tl);

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                if (i==1)
                    Tl_[cellI];// = T[cellI];//*0. + 300.;//T[cellI];

                // Calculate saturated density, specific humidity
                rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
                evsat_[cellI] = calc_evsat(Tl_[cellI]);
                wsat_[cellI] = 0.621945*(evsat_[cellI]/(p_-evsat_[cellI])); // ASHRAE 1, eq.23

                // Calculate transpiration rate]);
                E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(wsat_[cellI]-w[cellI])/(ra_[cellI]+rs_[cellI]);
                //E_[cellI] = 0.0; // No evapotranspiration

                // Calculate latent heat flux
                Ql_[cellI] = lambda_.value()*E_[cellI];

                // Calculate new leaf temperature
                new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Ql_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));

            }
        }

        // info
        Info << "    max leaf temp tl=" << gMax(new_Tl.internalField())
             << " K, iteration i="   << i << endl;

        // Check rel. L-infinity error
        maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
        maxRelError = maxError/gMax(mag(new_Tl.internalField()));

        // update leaf temp.
        forAll(Tl_, cellI)
            Tl_[cellI] = 0.5*Tl_[cellI]+0.5*new_Tl[cellI];

         // convergence check
         if (maxRelError < 1e-8)
             break;
    }
    Tl_.correctBoundaryConditions();

    // Iteration info
    Info << "Vegetation model:  Solving for Tl, Final residual = " << maxError
         << ", Final relative residual = " << maxRelError
         << ", No Iterations " << i << endl;

    Info << "temperature parameters: max Tl = " << gMax(Tl_)
         << ", min T = " << gMin(T) << ", max T = " << gMax(T) << endl;

    Info << "resistances: max rs = " << gMax(rs_)
         << ", max ra = " << gMax(ra_) << endl;

    // Final: Solve aerodynamc, stomatal resistance
    resistance(magU, T, w, Tl_);

    // Final: Update sensible and latent heat flux
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            // Calculate saturated density, specific humidity
            rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
            evsat_[cellI] = calc_evsat(Tl_[cellI]);
            wsat_[cellI] = 0.621945*(evsat_[cellI]/(p_-evsat_[cellI])); // ASHRAE 1, eq.23

            // Calculate transpiration rate
            E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(wsat_[cellI]-w[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side

            // Calculate latent heat flux
            Ql_[cellI] = lambda_.value()*E_[cellI];

            // Calculate sensible heat flux
            Qs_[cellI] = 2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]*(Tl_[cellI]-T[cellI])/ra_[cellI];
        }
    }
    rhosat_.correctBoundaryConditions();
    wsat_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    Ql_.correctBoundaryConditions();
    Qs_.correctBoundaryConditions();


}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> modifiedVegetationModel::ST()
{
    ST_ = Qs_/(rhoa_*cpa_);
    ST_.correctBoundaryConditions();
    return ST_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> modifiedVegetationModel::Su(volVectorField& U)
{
    Su_ = -Cf_*mag(U)*U;

    return fvm::SuSp(-Cf_*mag(U), U);
}

// return specific humidity source term
tmp<volScalarField> modifiedVegetationModel::Sw()
{
    Sw_ = E_/rhoa_;
    Sw_.correctBoundaryConditions();
    return Sw_;
}

// -----------------------------------------------------------------------------

bool modifiedVegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam
