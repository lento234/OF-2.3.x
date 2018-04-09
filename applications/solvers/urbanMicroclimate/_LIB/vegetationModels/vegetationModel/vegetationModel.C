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

################
Vegetation model implemented by L. Manickathan, Empa, February 2017
################      

\*---------------------------------------------------------------------------*/

#include "vegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(vegetationModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
vegetationModel::vegetationModel
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
    qsat_
    (
        IOobject
        (
            "qsat",
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
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double vegetationModel::calc_evsat(double& T)
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
double vegetationModel::calc_rhosat(double& T)
{
    return calc_evsat(T)/(461.5*T);
}

// solve radiation
void vegetationModel::radiation()
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
void vegetationModel::resistance(volScalarField& magU, volScalarField& T, volScalarField& q, volScalarField& Tl)
{
    const double p_ = 101325;

    // Calculate magnitude of velocity and bounding above Umin
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            //Aerodynamic resistance
            // ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);
            ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

            // Calculate vapor pressure of air
            //ev_[cellI] = q[cellI]*rhoa_.value()*T[cellI]*461.5;
            ev_[cellI] = p_*q[cellI]/(0.621945+q[cellI]);

            // Calculate sat. vapor pressure of air
            //evsat_[cellI] = calc_evsat(T[cellI]); // TODO bug
            evsat_[cellI] = calc_evsat(T[cellI]);

            // Vapor pressure deficit - kPa
            // VPD_[cellI] = (calc_evsat(T[cellI]) - (q[cellI]*rhoa_.value()*T[cellI]*461.5))/1000.0; // kPa
            //VPD_[cellI] = ev_[cellI] - evsat_[cellI];
            VPD_[cellI] = evsat_[cellI] - ev_[cellI];


            // Stomatal resistance - type 1
            // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
            //rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);
            // rs_[cellI] = rsMin_.value()*(31.0 + Rg_[cellI].component(2))*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rg_[cellI].component(2));


            // Stomatal resistance - type 2
            // rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2)); // type 2
            //if ((VPD_[cellI]/1000.0) < D0_.value())
            //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()));
            //else
            //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));
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
void vegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& q)
{
    // solve radiation within vegetation
    radiation();

    const double p_ = 101325;

    // Magnitude of velocity
    volScalarField magU("magU", mag(U));
    // Bounding velocity
    bound(magU, UMin_);

    // solve aerodynamic, stomatal resistance
    //resistance(U,T);
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
        resistance(magU, T, q, new_Tl);

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                if (i==1)
                    Tl_[cellI];// = T[cellI];//*0. + 300.;//T[cellI];

                // Calculate saturated density, specific humidity
                rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
                //qsat_[cellI]   = rhosat_[cellI]/rhoa_.value();
                evsat_[cellI] = calc_evsat(Tl_[cellI]);
                qsat_[cellI] = 0.621945*(evsat_[cellI]/(p_-evsat_[cellI])); // ASHRAE 1, eq.23

                // Calculate transpiration rate
                // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());
                // E_[cellI] = 2.0*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
                E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
                //E_[cellI] = 0.0; // No evapotranspiration

                // Calculate latent heat flux
                //Ql_[cellI] = lambda_.value()*E_[cellI];
                Ql_[cellI] = lambda_.value()*E_[cellI];

                // Calculate new leaf temperature
                //new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Ql_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));
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
        //Tl_.internalField() = new_Tl.internalField();
        //Tl_.boundaryField() = new_Tl.boundaryField();
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
    resistance(magU, T, q, Tl_);

    // Final: Update sensible and latent heat flux
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            // Calculate saturated density, specific humidity
            rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
            // qsat_[cellI] = rhosat_[cellI]/rhoa_.value();
            evsat_[cellI] = calc_evsat(Tl_[cellI]);
            qsat_[cellI] = 0.621945*(evsat_[cellI]/(p_-evsat_[cellI])); // ASHRAE 1, eq.23

            // Calculate transpiration rate
            // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());
            // E_[cellI] = 2.0*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
            E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side
            //E_[cellI] = 0.0; // no evapotranspiration
            // TODO: flag for no transpiration, one side, both side

            // Calculate latent heat flux
            Ql_[cellI] = lambda_.value()*E_[cellI];

            // Calculate sensible heat flux
            Qs_[cellI] = 2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]*(Tl_[cellI]-T[cellI])/ra_[cellI];
        }
    }
    rhosat_.correctBoundaryConditions();
    qsat_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    Ql_.correctBoundaryConditions();
    Qs_.correctBoundaryConditions();

    // Iteration info
    // Info << "              Vegetation model:  max. Rn = " << max(mag(Rn_))
    //      << "; max. Ql = " << max(mag(Ql_))
    //      << "; max. Qs = " << max(mag(Qs_))
    //      << "; error: max. Esum = " << max(mag(Rn_.internalField() - Qs_.internalField()- Ql_.internalField())) << endl;

}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> vegetationModel::Sh()
{
    // forAll(LAD_, cellI)
    //     if (LAD_[cellI] > 10*SMALL)
    //         Sh_[cellI] = Qs_[cellI]/(rhoa_.value()*cpa_.value());
    // Sh_.correctBoundaryConditions();
    // return Sh_;
    Sh_ = Qs_/(rhoa_*cpa_);
    Sh_.correctBoundaryConditions();
    return Sh_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> vegetationModel::Su(volVectorField& U)
{
    // forAll(Su_, cellI)
    //     Su_[cellI] = -Cf_[cellI]*mag(U[cellI])*U[cellI]; // Cf = LAD*Cd
    // //Su_[cellI] = -0.5*Cf_[cellI]*mag(U[cellI])*U[cellI];
    // Su_.correctBoundaryConditions();
    // return Su_;
    return fvm::SuSp(-Cf_*mag(U), U);
}

// return specific humidity source term
tmp<volScalarField> vegetationModel::Sq()
{
    // forAll (LAD_, cellI)
    //     if (LAD_[cellI] > 10*SMALL)
    //         Sq_[cellI] = E_[cellI];
    // Sq_.correctBoundaryConditions();
    // return Sq_;
    Sq_ = E_/rhoa_;
    Sq_.correctBoundaryConditions();
    return Sq_;
}

// -----------------------------------------------------------------------------

bool vegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam



// tmp<volScalarField> vegetationModel::Sh()
// {
//     forAll(LAD_, cellI)
//         if (LAD_[cellI] > 10*SMALL)
//             Sh_[cellI] = Qs_[cellI]/(rhoa_.value()*cpa_/value());
//
//     return Sh_;
//     // return fvm::SuSp((Qs(U,T)/(rhoa_*cpa_))/T,T);
// }
// void vegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& q)
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


// void vegetationModel::testing(volVectorField& U, volScalarField& T, volScalarField& q)
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
// double vegetationModel::calc_qsat(double& rhosat)
// {
//     return rhosat/rhoa_.value();
// }
//
// // calc tranpiration rate
// double vegetationModel::calc_E(doublvolVectorField& U, volScalarField& T, volScalarField& q)
// {
//     return LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_)
//     volScalarField Tleaf("Tleaf", Tl_);
//     return tmp<volScalarField>
//     (
//         new volScalarField("E", LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_))
//     );
// }

// return tranpiration rate
// tmp<volScalarField> vegetationModel::E(volVectorField& U, volScalarField& T, volScalarField& q) const
// {
//     volScalarField Tleaf("Tleaf", Tl_);
//     return tmp<volScalarField>
//     (
//         new volScalarField("E", LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_))
//     );
// }

// return latent heat
// tmp<volScalarField> vegetationModel::Ql(volVectorField& U, volScalarField& T, volScalarField& q) const
// {
//     return tmp<volScalarField>
//     (
//         new volScalarField("Ql", lambda_*E(U,T,q))
//     );
// }

// return sensible heat
// tmp<volScalarField> vegetationModel::Qs(volVectorField& U, volScalarField& T) const
// {
//     return tmp<volScalarField>
//     (
//         new volScalarField("Qs", 2.0*rhoa_*cpa_*LAD_*(Tl_-T)/ra(U))
//     );
// }



// tmp<volScalarField> vegetationModel::rhosat(volScalarField& T) const
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
// tmp<volScalarField> vegetationModel::Rg() const
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
// tmp<volVectorField> vegetationModel::Rg()
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

// tmp<volScalarField> vegetationModel::Rg()
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
// tmp<volScalarField> vegetationModel::Rn()
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
// tmp<volScalarField> vegetationModel::qsat(volScalarField& Tl) const
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

// void vegetationModel::testing(volVectorField& U, volScalarField& T, volScalarField& q)
// {
//     update_Tl(U, T, q);
// }

// return temp source term
// tmp<fvScalarMatrix> vegetationModel::Sh(volVectorField& U, volScalarField& T) const
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
