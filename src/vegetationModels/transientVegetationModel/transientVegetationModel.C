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

#include "transientVegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(transientVegetationModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transientVegetationModel::transientVegetationModel
(
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& Yv
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
        TimeDataEntry<scalar>
        (
            runTime_.time(),
            "Rg0",
            vegetationProperties_
        )
    ),
    Rl0_
    (
        TimeDataEntry<scalar>
        (
            runTime_.time(),
            "Rl0",
            vegetationProperties_
        )
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
    pv_
    (
        IOobject
        (
            "pv",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),
    pvsat_
    (
        IOobject
        (
            "pvsat",
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
    Yvsat_
    (
        IOobject
        (
            "Yvsat",
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
    SYv_
    (
        IOobject
        (
            "SYv",
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

        // Extract value
        a1 = a1_.value();
        a2 = a2_.value();
        a3 = a3_.value();
        cpa = cpa_.value();
        C = C_.value();
        D0 = D0_.value();
        H = H_.value();
        kc = kc_.value();
        l = l_.value();
        rhoa = rhoa_.value();
        rsMin = rsMin_.value();
        lambda = lambda_.value();

        Info << " Defined custom vegetation model" << endl;
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// calc saturated water vapor pressure
double transientVegetationModel::calc_pvsat(double& T)
{
    return exp( 77.3450 + 0.0057*T - 7235.0/T) / pow(T, 8.2);
}

// calc saturated density of water vapour
double transientVegetationModel::calc_rhosat(double& T)
{
    //return 0.0022 * exp( 77.3450 + 0.0057*T - 7235.0/T) / pow(T, 9.2);
    return calc_pvsat(T)/(461.5*T);
}

// solve radiation
void transientVegetationModel::radiation()
{
    // Fetch time dependent radiation parameters
    Rg0 = Rg0_.value(runTime_.value());
    Rl0 = Rl0_.value(runTime_.value());

    // empirical global radiation = f(height) inside and below vegetation
    forAll(Rg_, cellI)
        Rg_[cellI] = vector(0,0, Rg0*exp(-kc*LAI_[cellI]));
    // Rg_[cellI] = vector(0,0, Rg0_.value(runTime_.value())*exp(-kc_.value()*LAI_[cellI]));
    // Rg_[cellI] = vector(0,0, Rg0*exp(-kc*LAI_[cellI]));
    //  Rg_[cellI] = vector(0,0, Rg0_.value()*exp(-kc_.value()*LAI_[cellI]));
    //  Rg_[cellI] = vector(0,0, Rg0_.value()*exp(-kc_.value()*LAI_[cellI]*
    //                             (H_.value() - mesh_.C()[cellI].component(2))/H_.value()));
    Rg_.correctBoundaryConditions();

    volTensorField gradRg = fvc::grad(Rg_);

    // radiation density inside vegetation
    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Rn_[cellI] = gradRg[cellI].component(8) + 0.04*(Rl0/H);
    // Rn_[cellI] = gradRg[cellI].component(8) + 0.04*(Rl0_.value(runTime_.value())/H_.value());
    // Rn_[cellI] = gradRg[cellI].component(8) + 0.04*(Rl0/H_.value());
    // Rn_[cellI] = gradRg[cellI].component(8) + 0.04*(Rl0_.value()/H_.value()); //gradRg[cellI] && tensor(0,0,0,0,0,0,0,0,1);
    // Rn_[cellI] = gradRg[cellI].component(8);
    Rn_.correctBoundaryConditions();

}

// solve aerodynamic resistance
void transientVegetationModel::resistance(volScalarField& magU, volScalarField& T, volScalarField& Yv)
{

    // Calculate magnitude of velocity and bounding above Umin
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            //Aerodynamic resistance
            ra_[cellI] = C*pow(l/magU[cellI], 0.5);
            // ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);
            // ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

            // Calculate vapor pressure
            pv_[cellI] = Yv[cellI]*rhoa*T[cellI]*461.5; // TODO: not correct
            // TODO: not correct
            // q = Yv[cellI]/(1-Yv[cellI]);
            // pv_[cellI] = q*rhoa*T[cellI]*461.5;
            // ev_[cellI] = q[cellI]*rhoa_.value()*T[cellI]*461.5;

            // Calculate sat. vapor pressure
            pvsat_[cellI] = calc_pvsat(T[cellI]); // TODO: functiion of leaf temperature, not air.

            // Vapor pressure deficit - kPa
            // VPD_[cellI] = (calc_evsat(T[cellI]) - (q[cellI]*rhoa_.value()*T[cellI]*461.5))/1000.0; // kPa
            // VPD_[cellI] = pv_[cellI] - pvsat_[cellI];
            VPD_[cellI] = (pv_[cellI] - pvsat_[cellI]);

            // Stomatal resistance - type 1
            // rs_[cellI] = rsMin*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
            //rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
            //rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);


            // Stomatal resistance - type 2
            rs_[cellI] = rsMin*((a1 + Rg0)/(a2 + Rg0))*(1.0 + a3*pow(VPD_[cellI]/1000.0-D0,2));
            // rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2)); // type 2

            // Stomatal resistance - type 2b
            if ((VPD_[cellI]/1000.0) < D0_.value())
                rs_[cellI] = rsMin*((a1 + Rg0)/(a2 + Rg0));
            else
                rs_[cellI] = rsMin*((a1 + Rg0)/(a2 + Rg0))*(1.0 + a3*pow(VPD_[cellI]/1000.0-D0,2));
        }
    }
    pv_.correctBoundaryConditions();
    pvsat_.correctBoundaryConditions();
    VPD_.correctBoundaryConditions();
    ra_.correctBoundaryConditions();
    rs_.correctBoundaryConditions();
}

// solve vegetation model
void transientVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& Yv)
{
    // solve radiation within vegetation
    radiation();

    // Magnitude of velocity
    volScalarField magU("magU", mag(U));
    // Bounding velocity
    bound(magU, UMin_);

    // solve aerodynamic, stomatal resistance
    //resistance(U,T);
    volScalarField new_Tl("new_Tl", Tl_);

    // solve leaf temperature, iteratively.
    int maxIter = 100;
    for (int i=1; i<=maxIter; i++)
    {
        // Solve aerodynamc, stomatal resistance
        resistance(magU, new_Tl, Yv);

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                // if (i==0)
                    // Tl_[cellI] = T[cellI];

                // Calculate saturated density, specific humidity
                rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
                // qsat_[cellI]   = rhosat_[cellI]/rhoa_.value();
                Yvsat_[cellI]   = rhosat_[cellI]/rhoa; // TODO: Yvsat_ = qsat/(1+qsat)
                // qsat_[cellI] = rhosat_[cellI]/rhoa_.value();
                // Yvsat_[cellI] = qsat_[cellI]/(1+qsat_[cellI]);
                // TODO: Yvsat_ = qsat/(1+qsat)

                // Calculate transpiration rate
                // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());
                // E_[cellI] = 2.0*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
                // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
                E_[cellI] = LAD_[cellI]*rhoa*(Yvsat_[cellI]-Yv[cellI])/(ra_[cellI]+rs_[cellI]);

                // Calculate latent heat flux
                // Ql_[cellI] = lambda_.value()*E_[cellI];
                Ql_[cellI] = lambda*E_[cellI];

                // Calculate new leaf temperature
                // new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Ql_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));
                new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Ql_[cellI])*(ra_[cellI]/(2.0*rhoa*cpa*LAD_[cellI]));
            }
        }
        new_Tl.correctBoundaryConditions();

        // Check rel. L-infinity error
        scalar maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
        scalar maxRelError = maxError/gMax(mag(new_Tl.internalField()));

        // Iteration info
        Info << "      Vegetation model:  Solving for Tl. Iteration " << i
             << "; Rg0 = " << Rg0_.value(runTime_.value())
             << "; Rl0 = " << Rl0_.value(runTime_.value())
             << "; max. error = " << maxError
             << "; max. rel. error = " << maxRelError << endl;;


        // update leaf temp.
        Tl_.internalField() = new_Tl.internalField();
        Tl_.boundaryField() = new_Tl.boundaryField();

         // convergence check
         if (maxRelError < 1e-8)
         {
             Info << "      Vegetation model:  CONVERGED" << endl;
             break;
         }

         if (i == maxIter)
             Info << "      Vegetation model:  >>>>>>>>>>>>>>> N O T  C O N V E R G E D  !! <<<<<<<<<<<<<" << endl;
    }

    // Final: Solve aerodynamc, stomatal resistance
    resistance(magU, new_Tl, Yv);

    // Final: Update sensible and latent heat flux
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            // Calculate saturated density, specific humidity
            rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
            // qsat_[cellI] = rhosat_[cellI]/rhoa_.value();
            Yvsat_[cellI] = rhosat_[cellI]/rhoa; // TODO: not correct
            // qsat_[cellI] = rhosat_[cellI]/rhoa_.value();
            // Yvsat_[cellI] = qsat_[cellI]/(1+qsat_[cellI]);
            // TODO: Yvsat_ = qsat/(1+qsat)

            // Calculate transpiration rate
            // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());
            // E_[cellI] = 2.0*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
            // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side
            E_[cellI] = LAD_[cellI]*rhoa*(Yvsat_[cellI]-Yv[cellI])/(ra_[cellI]+rs_[cellI]); // TODO: implement switch for double or single side

            // Calculate latent heat flux
            // Ql_[cellI] = lambda_.value()*E_[cellI];
            Ql_[cellI] = lambda*E_[cellI];

            // Calculate sensible heat flux
            // Qs_[cellI] = 2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]*(Tl_[cellI]-T[cellI])/ra_[cellI];
            Qs_[cellI] = 2.0*rhoa*cpa*LAD_[cellI]*(Tl_[cellI]-T[cellI])/ra_[cellI];
        }
    }
    rhosat_.correctBoundaryConditions();
    Yvsat_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    Ql_.correctBoundaryConditions();
    Qs_.correctBoundaryConditions();
    //
    // scalar curTime = runTime_.value();
    // Info << "runTime_" << curTime
    //      << ", testVar: " << testVar.value(curTime) << endl;
    // Rg0_.value(runTime_.value())
    // Iteration info
    // Info << "              Vegetation model:  max. Rn = " << max(mag(Rn_))
    //      << "; max. Ql = " << max(mag(Ql_))
    //      << "; max. Qs = " << max(mag(Qs_))
    //      << "; error: max. Esum = " << max(mag(Rn_.internalField() - Qs_.internalField()- Ql_.internalField())) << endl;

}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> transientVegetationModel::Sh()
{
    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Sh_[cellI] = Qs_[cellI]/(rhoa*cpa);
    // Sh_[cellI] = Qs_[cellI]/(rhoa_.value()*cpa_.value());
    Sh_.correctBoundaryConditions();
    return Sh_;
}

// solve & return momentum source term (explicit)
tmp<volVectorField> transientVegetationModel::Su(volVectorField& U)
{
    forAll(Su_, cellI)
        Su_[cellI] = -Cf_[cellI]*mag(U[cellI])*U[cellI]; // Cf = LAD*Cd
    //Su_[cellI] = -0.5*Cf_[cellI]*mag(U[cellI])*U[cellI];
    Su_.correctBoundaryConditions();
    return Su_;
}

// return specific humidity source term
tmp<volScalarField> transientVegetationModel::SYv()
{
    forAll (LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            SYv_[cellI] = E_[cellI];
    SYv_.correctBoundaryConditions();
    return SYv_;
}

// -----------------------------------------------------------------------------

bool transientVegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam



// tmp<volScalarField> transientVegetationModel::Sh()
// {
//     forAll(LAD_, cellI)
//         if (LAD_[cellI] > 10*SMALL)
//             Sh_[cellI] = Qs_[cellI]/(rhoa_.value()*cpa_/value());
//
//     return Sh_;
//     // return fvm::SuSp((Qs(U,T)/(rhoa_*cpa_))/T,T);
// }
// void transientVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& q)
// {
//     // net radiation
//     volScalarField tmp_Rn("tmp_Rn", Rn(U));
//
//     // vegetation resistance
//     volScalarField tmp_ra("tmp_ra", ra(U));
//
//     // bounding LAD
//     volScalarField tmp_LAD("tmp_LAD", LAD_);
//     dimensionedScalar LADMin("LADMin", dimensionSet(0,-1,0,0,0,0,0), SMALL);
//     bound(tmp_LAD, LADMin);
//
//     // Initial leaf temperature
//     // forAll(tmp_LAD, cellI)
//     // {
//     //     if (tmp_LAD[cellI] > SMALL)
//     //         Tl[cellI] = T[cellI];
//     // }
//     Tl_.internalField() = T.internalField();
//
//     // Info << "WIP <<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
//     // solve for leaf temperature
//     int i = 0; int maxIter = 10;
//     while (i < maxIter)
//     {
//         // Calculate latent heat flux.
//         volScalarField tmp_Ql("tmp_Ql", Ql(U, T, q));
//
//         //Calculate new leaf temperature
//         volScalarField new_Tl("new_Tl", T + (tmp_Rn - tmp_Ql)*(tmp_ra/(2.0*rhoa_*cpa_*tmp_LAD)) );
//
//         // Determine errors
//         scalar maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
//         scalar maxRelError = maxError/gMax(mag(new_Tl.internalField()));
//         // Iteration info
//         Info << "Calculating leaf temp. Iteration i : " << i++
//              << ", max. error: " << maxError
//              << ", max. rel. error: " << maxRelError << endl;;
//
//         // update leaf temp.
//         Tl_.internalField() = new_Tl.internalField();
//         // forAll(tmp_LAD, cellI)
//         // {
//         //     if (tmp_LAD[cellI] > 100*SMALL)
//         //         Tl_.internalField()[cellI] = new_Tl.internalField()[cellI];
//         // }
//
//         // convergence check
//         if (maxRelError < 1e-8)
//         {
//             Info << "Leaf temperature calculation CONVERGED !!"
//                  << ", max. rel. error: " << maxRelError << endl;
//             break;
//         }
//
//     }
//
//     // Cleanup Tl at region without leaf
//     forAll(tmp_LAD, cellI)
//     {
//         if (tmp_LAD[cellI] <= 10*SMALL)
//             Tl_[cellI] = gMin(Tl_);
//     }
//
//     if (i == maxIter)
//     {
//         Info << "WARNING!! Leaf temperature not converged !!" << endl;
//     }
//
//     Info << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Runtime time: " << runTime_.timeName() << endl;
// }


// void transientVegetationModel::testing(volVectorField& U, volScalarField& T, volScalarField& q)
// {
//     // volVectorField tmpRg("Rg", Rg());
//     // volScalarField tmpRn("Rn", Rn());
//
//     // radiation();
//
//     // Info << ">>>>>>>> max Rg: " << gMax(tmpRg) << endl;
//     // Info << ">>>>>>>> max Rn: " << gMax(tmpRn) << endl;
//     Info << "Hi" << endl;
// }

// // calc saturated specific humidity of water vapour
// double transientVegetationModel::calc_qsat(double& rhosat)
// {
//     return rhosat/rhoa_.value();
// }
//
// // calc tranpiration rate
// double transientVegetationModel::calc_E(doublvolVectorField& U, volScalarField& T, volScalarField& q)
// {
//     return LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_)
//     volScalarField Tleaf("Tleaf", Tl_);
//     return tmp<volScalarField>
//     (
//         new volScalarField("E", LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_))
//     );
// }

// return tranpiration rate
// tmp<volScalarField> transientVegetationModel::E(volVectorField& U, volScalarField& T, volScalarField& q) const
// {
//     volScalarField Tleaf("Tleaf", Tl_);
//     return tmp<volScalarField>
//     (
//         new volScalarField("E", LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_))
//     );
// }

// return latent heat
// tmp<volScalarField> transientVegetationModel::Ql(volVectorField& U, volScalarField& T, volScalarField& q) const
// {
//     return tmp<volScalarField>
//     (
//         new volScalarField("Ql", lambda_*E(U,T,q))
//     );
// }

// return sensible heat
// tmp<volScalarField> transientVegetationModel::Qs(volVectorField& U, volScalarField& T) const
// {
//     return tmp<volScalarField>
//     (
//         new volScalarField("Qs", 2.0*rhoa_*cpa_*LAD_*(Tl_-T)/ra(U))
//     );
// }



// tmp<volScalarField> transientVegetationModel::rhosat(volScalarField& T) const
// {
//     volScalarField rhosat
//     (
//         IOobject
//         (
//             "rhosat",
//             T.time().timeName(),
//             T.mesh(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         T.mesh(),
//         dimensionedScalar("0", dimDensity, 0.0)
//     );
//
//     // Calculate saturated water vapour density
//     rhosat.internalField() = 0.0022* exp(  77.3450 + 0.0057*T.internalField()
//                                          - 7235.0/T.internalField()) / pow(T.internalField(), 9.2);
//
//     return tmp<volScalarField>
//     (
//         new volScalarField("rhosat", rhosat)
//     );
// }

// return empirical radiation
// tmp<volScalarField> transientVegetationModel::Rg() const
// {
//     volScalarField Rg
//     (
//         IOobject
//         (
//             "Rg",
//             runTime_.timeName(),
//             mesh_,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh_,
//         dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
//     );
//
//     forAll(Rg, cellI)
//     {
//         if (LAD_ > 10*SMALL)
//             Rg[cellI] = Rg0_*LAI_*exp(kc_*(H_ - U.mesh().C().component(2))/H_)
//     }
//
//     return tmp<volScalarField>
//     (
//         new volScalarField("Rg", Rg0_*LAI_*exp(kc_*(H_ - U.mesh().C().component(2))/H_))
//     );
// }

// return empirical radiation
// tmp<volVectorField> transientVegetationModel::Rg()
// {
//     forAll(LAD_, cellI)
//     {
//         if (LAD_[cellI] > 10*SMALL)
//             Rg_[cellI] = vector(0,0, Rg0_.value()*LAI_[cellI] * exp(kc_.value() *
//                                      (H_.value() - mesh_.C()[cellI].component(2))/H_.value()));
//     }
//
//     return Rg_;
// }



// Info << "Rg: " << (fvc::grad(Rg_[cellI] & vector(0,0,1))) << endl;
// Rn_[cellI] = fvc::grad(Rg_[cellI]) && tensor(0,0,0,0,0,0,0,0,1);

// tmp<volScalarField> transientVegetationModel::Rg()
// {
//     forAll(LAD_, cellI)
//     {
//         if (LAD_[cellI] > 10*SMALL)
//             Rg_[cellI] = Rg0_.value()*LAI_[cellI] * exp(kc_.value() *
//                             (H_.value() - mesh_.C()[cellI].component(2))/H_.value());
//     }
//
//     return Rg_;
// }

// return net radiation density in volume
// tmp<volScalarField> transientVegetationModel::Rn()
// {
//     // double tmpRg;
//     // forAll(LAD_, cellI)
//     // {
//     //     if (LAD_[cellI] > 10*SMALL)
//     //     {
//     //         // Radiation
//     //         Rg_[cellI] = Rg0_.value()*LAI_[cellI] * exp(kc_.value() *
//     //                         (H_.value() - mesh_.C()[cellI].component(2))/H_.value());
//     //     }
//     //         Rn_[cellI] = fvc::grad(Rg()) & vector(0,0,1);
//     // }
//     Rn_.internalField() = fvc::grad(Rg()) & vector(0,0,1);
//     // return tmp<volScalarField>
//     // (
//         // new volScalarField("Rn", fvc::grad(Rg()) & vector(0,0,1))
//     // );
//     return Rn_;
// }

// return saturated specific humidity of water vapour
// tmp<volScalarField> transientVegetationModel::qsat(volScalarField& Tl) const
// {
//     volScalarField qsat
//     (
//         IOobject
//         (
//             "qsat",
//             Tl.time().timeName(),
//             Tl.mesh(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         Tl.mesh(),
//         dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0)
//     );
//
//     // Calculate satura
//     qsat.internalField() = rhosat(Tl)/rhoa_;
//
//     return tmp<volScalarField>
//     (
//         new volScalarField("qsat", qsat)
//     );
// }

// ************************************************************************* //

// void transientVegetationModel::testing(volVectorField& U, volScalarField& T, volScalarField& q)
// {
//     update_Tl(U, T, q);
// }

// return temp source term
// tmp<fvScalarMatrix> transientVegetationModel::Sh(volVectorField& U, volScalarField& T) const
// {
//     // Calculate magnitude of velocity and bounding above Umin
//     volScalarField magU("magU", mag(U));
//     bound(magU, UMin_);
//
//     // calculate aerodynamic resitance
//     volScalarField ra("ra", C_*pow(l_/magU, 0.5));
//
//     return
//     (
//         fvm::SuSp((2*a_*(Tl_-T)/ra)/T,T)
//     );
// }

//
//
// Info << "newTl: " << max(newTl) << endl;
// Info << "Info: " << max(tmp_ra) << endl;
// Info << "Info: " << max(rhosat(T)) << endl;
// Info << "Info: " << max(qsat(T)) << endl;

// volScalarField rhosat
// (
//     IOobject
//     (
//         "rhosat",
//         T.time().timeName(),
//         T.mesh(),
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     U.mesh(),
//     dimensionedScalar("0", dimDensity, 0.0)
// );
//
// // Calculate satura
// rhosat.internalField() = exp(  77.3450 + 0.0057*T.internalField()
//                              - 7235.0/T.internalField()) / pow(T.internalField(), 9.2);

// double tempT;
// forAll(rhosat, cellI)
// {
//     tempT = T[cellI];
//     rhosat[cellI] = exp(77.3450 + 0.0057*tempT - 7235.0/tempT)/pow(tempT,9.2);
// }
// dimensionedScalar("0", dimensionSet(0, 2, -1 , 0, 0), 0.0)
//     // // solve for new Tl
//     // int i = 0;
//     // do
//     // {
//     //     volScalarField oldTl("oldTl", Tl_);
//     //     volScalarField oldQl("oldQl", Ql(U, T, q, oldTl));
//     //
//     //     volScalarField newTl("newTl", T + (Rnet - oldQl)*())
//     //
//     //     i++;
//     // } while( i < 5);
//
//     double normT = 0;
//     forAll(T, cellI)
//     {
//         Info << pow(T[cellI],2) * U.mesh().V()[cellI] << endl;
//         normT += pow(T[cellI],2) * U.mesh().V()[cellI];
//     }
//     Info << "Info: " << normT << endl;
// Info << "Info: " << max(tmp_Rn) << endl;
// Info << "Info: " << max(tmp_ra) << endl;
// Info << "Info: " << max(rhosat(T)) << endl;
// Info << "Info: " << max(qsat(T)) << endl;





    // Info << "Hi" << endl;

    // Info << rhosat(T[0]) << endl;
    // for (int i=0; i<20; i++)
    // {
    //     forAll(LAD_, cellI)
    //     {
    //         if (LAD_ > 10*SMALL)
    //         {
    //             // Initial leaf temperature
    //             if i==0:
    //                 Tl_[cellI] = T[cellI];
    //
    //             // Calculate latent heat flux.
    //             // Ql_[cellI] = calcQl(U[cellI]);
    //             Info <<
    //
    //             // Calculate new leaf temperature
    //
    //             // Check for convergence
    //
    //         }
    //     }
    // }


    // Info << tempData << endl;
    // net radiation
    // volScalarField tmp_Rn("tmp_Rn", Rn(U));

    // // vegetation resistance
    // volScalarField tmp_ra("tmp_ra", ra(U));
    //
    // // bounding LAD
    // volScalarField tmp_LAD("tmp_LAD", LAD_);
    // dimensionedScalar LADMin("LADMin", dimensionSet(0,-1,0,0,0,0,0), SMALL);
    // bound(tmp_LAD, LADMin);
    //
    // // Initial leaf temperature
    // // forAll(tmp_LAD, cellI)
    // // {
    // //     if (tmp_LAD[cellI] > SMALL)
    // //         Tl[cellI] = T[cellI];
    // // }
    // Tl_.internalField() = T.internalField();
    //
    // // Info << "WIP <<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    // // solve for leaf temperature
    // int i = 0; int maxIter = 10;
    // while (i < maxIter)
    // {
    //     // Calculate latent heat flux.
    //     volScalarField tmp_Ql("tmp_Ql", Ql(U, T, q));
    //
    //     //Calculate new leaf temperature
    //     volScalarField new_Tl("new_Tl", T + (tmp_Rn - tmp_Ql)*(tmp_ra/(2.0*rhoa_*cpa_*tmp_LAD)) );
    //
    //     // Determine errors
    //     scalar maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
    //     scalar maxRelError = maxError/gMax(mag(new_Tl.internalField()));
    //     // Iteration info
    //     Info << "Calculating leaf temp. Iteration i : " << i++
    //          << ", max. error: " << maxError
    //          << ", max. rel. error: " << maxRelError << endl;;
    //
    //     // update leaf temp.
    //     Tl_.internalField() = new_Tl.internalField();
    //     // forAll(tmp_LAD, cellI)
    //     // {
    //     //     if (tmp_LAD[cellI] > 100*SMALL)
    //     //         Tl_.internalField()[cellI] = new_Tl.internalField()[cellI];
    //     // }
    //
    //     // convergence check
    //     if (maxRelError < 1e-8)
    //     {
    //         Info << "Leaf temperature calculation CONVERGED !!"
    //              << ", max. rel. error: " << maxRelError << endl;
    //         break;
    //     }
    //
    // }
    //
    // // Cleanup Tl at region without leaf
    // forAll(tmp_LAD, cellI)
    // {
    //     if (tmp_LAD[cellI] <= 10*SMALL)
    //         Tl_[cellI] = gMin(Tl_);
    // }
    //
    // if (i == maxIter)
    // {
    //     Info << "WARNING!! Leaf temperature not converged !!" << endl;
    // }
    //
    // Info << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Runtime time: " << runTime_.timeName() << endl;
