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
    const volScalarField& LAD,
    const volScalarField& LAI,
    const volScalarField& T
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
    C_
    (
        vegetationProperties_.lookup("C")
    ),
    Cdf_
    (
        vegetationProperties_.lookup("Cdf")
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
    rhoa_
    (
        vegetationProperties_.lookup("rhoa")
    ),
    rs_
    (
        vegetationProperties_.lookup("rs")
    ),
    TlMin_("TlMin", dimTemperature, SMALL),
    UMin_("UMin", dimVelocity, SMALL),
    lambda_
    (
        vegetationProperties_.lookup("lambda")
    ),
    LAD_(LAD),
    LAI_(LAI),
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
    )
    {
        // Bounding parameters
        bound(Tl_, TlMin_);

        Info << "Defined custom vegetation model" << endl;
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// calc saturated density of water vapour
double vegetationModel::calc_rhosat(double& T)
{
    return 0.0022 * exp( 77.3450 + 0.0057*T - 7235.0/T) / pow(T, 9.2);
}

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

// solve radiation
void vegetationModel::radiation()
{
    // empirical global radiation = f(height) inside vegetation
    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Rg_[cellI] = vector(0,0, Rg0_.value()*LAI_[cellI] * exp(kc_.value() *
                                        (H_.value() - mesh_.C()[cellI].component(2))/H_.value()));

    volTensorField gradRg = fvc::grad(Rg_);

    // radiation density inside vegetation
    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Rn_[cellI] = gradRg[cellI] && tensor(0,0,0,0,0,0,0,0,1);

}

// solve aerodynamic resistance
void vegetationModel::resistance(volVectorField& U)
{
    // Calculate magnitude of velocity and bounding above Umin
    volScalarField magU("magU", mag(U));
    bound(magU, UMin_);

    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

}




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

// -----------------------------------------------------------------------------

// return energy source term
// tmp<volScalarField> vegetationModel::Sh()
// {
//     forAll(LAD_, cellI)
//         if (LAD_[cellI] > 10*SMALL)
//             Sh_[cellI] = Qs_[cellI]/(rhoa_.value()*cpa_/value());
//
//     return Sh_;
//     // return fvm::SuSp((Qs(U,T)/(rhoa_*cpa_))/T,T);
// }

// solve & return momentum source term (explicit)
tmp<volVectorField> vegetationModel::Su(volVectorField& U)
{
    forAll(LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Su_[cellI] = -0.5*Cdf_.value()*LAD_[cellI]*mag(U[cellI])*U[cellI];

    return Su_;
}

// return specific humidity source term
tmp<volScalarField> vegetationModel::Sq()
{
    forAll (LAD_, cellI)
        if (LAD_[cellI] > 10*SMALL)
            Sq_[cellI] = E_[cellI];

    return Sq_;
}

// -----------------------------------------------------------------------------

// solve vegetation model
void vegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& q)
{
    // solve radiation within vegetation
    radiation();

    // solve aerodynamic, stomatal resistance
    resistance(U);

    int maxIter = 10;
    for (int i=0; i<maxIter; i++)
    {
        volScalarField new_Tl("new_Tl", Tl_);

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                if (i==0)
                    Tl_[cellI] = T[cellI];

                // Calculate saturated density, specific humidity
                rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
                qsat_[cellI]   = rhosat_[cellI]/rhoa_.value();

                // Calculate transpiration rate
                E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());

                // Calculate latent heat flux
                Ql_[cellI] = lambda_.value()*E_[cellI];

                // Calculate new leaf temperature
                new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Ql_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));
            }
        }

        // Check rel. L-infinity error
        scalar maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
        scalar maxRelError = maxError/gMax(mag(new_Tl.internalField()));

        // Iteration info
        Info << "       Calculating leaf temp. Iteration i : " << i
             << ", max. error: " << maxError
             << ", max. rel. error: " << maxRelError << endl;;


        // update leaf temp.
        Tl_.internalField() = new_Tl.internalField();

         // convergence check
         if (maxRelError < 1e-8)
         {
             Info << "          Leaf temperature calculation CONVERGED !!"
                  << ", max. rel. error: " << maxRelError << endl;
             break;
         }

         if (i == maxIter-1)
            Info << "           WARNING!! Leaf temperature NOT converged !!" << endl;
    }

}


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


bool vegetationModel::read()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam

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
