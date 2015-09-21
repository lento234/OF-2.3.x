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
    UMin_("UMin", dimVelocity, SMALL),
    lambda_
    (
        vegetationProperties_.lookup("lambda")
    ),
    LAD_(LAD),
    LAI_(LAI),
    Tl_(
        IOobject
        (
            "Tl",
            T.time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedScalar("0", dimensionSet(0,0,0,1,0,0,0), 0.0)
    )
    {
        Info << "Defined custom vegetation model" << endl;

        // Bounding Leaf temperature
        dimensionedScalar TlMin("TlMin", dimensionSet(0,0,0,1,0,0,0), SMALL);
        bound(Tl_, TlMin);
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// return tranpiration rate
tmp<volScalarField> vegetationModel::E(volVectorField& U, volScalarField& T, volScalarField& q) const
{
    volScalarField Tleaf("Tleaf", Tl_);
    return tmp<volScalarField>
    (
        new volScalarField("E", LAD_*rhoa_*(qsat(Tleaf)-q)/(ra(U)+rs_))
    );
}

// return latent heat
tmp<volScalarField> vegetationModel::Ql(volVectorField& U, volScalarField& T, volScalarField& q) const
{
    return tmp<volScalarField>
    (
        new volScalarField("Ql", lambda_*E(U,T,q))
    );
}

// return sensible heat
tmp<volScalarField> vegetationModel::Qs(volVectorField& U, volScalarField& T) const
{
    return tmp<volScalarField>
    (
        new volScalarField("Qs", 2.0*rhoa_*cpa_*LAD_*(Tl_-T)/ra(U))
    );
}

// return aerodynamic resistance
tmp<volScalarField> vegetationModel::ra(volVectorField& U) const
{
    // Calculate magnitude of velocity and bounding above Umin
    volScalarField magU("magU", mag(U));
    bound(magU, UMin_);

    return tmp<volScalarField>
    (
        new volScalarField("ra", C_*pow(l_/magU, 0.5))
    );

}

// return saturated density of water vapour
tmp<volScalarField> vegetationModel::rhosat(volScalarField& T) const
{
    volScalarField rhosat
    (
        IOobject
        (
            "rhosat",
            T.time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T.mesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    );

    // Calculate saturated water vapour density
    rhosat.internalField() = 0.0022* exp(  77.3450 + 0.0057*T.internalField()
                                         - 7235.0/T.internalField()) / pow(T.internalField(), 9.2);

    return tmp<volScalarField>
    (
        new volScalarField("rhosat", rhosat)
    );
}

// return empirical radiation
tmp<volScalarField> vegetationModel::Rg(volVectorField& U) const
{
    return tmp<volScalarField>
    (
        new volScalarField("Rg", Rg0_*LAI_*exp(kc_*(H_ - U.mesh().C().component(2))/H_))
    );
}

// return net radiation density in volume
tmp<volScalarField> vegetationModel::Rn(volVectorField& U) const
{
    return tmp<volScalarField>
    (
        new volScalarField("Rn", fvc::grad(Rg(U)) & vector(0,0,1))
    );
}

// return saturated specific humidity of water vapour
tmp<volScalarField> vegetationModel::qsat(volScalarField& Tl) const
{
    volScalarField qsat
    (
        IOobject
        (
            "qsat",
            Tl.time().timeName(),
            Tl.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Tl.mesh(),
        dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    // Calculate satura
    qsat.internalField() = rhosat(Tl)/rhoa_;

    return tmp<volScalarField>
    (
        new volScalarField("qsat", qsat)
    );
}

// -----------------------------------------------------------------------------

// return energy source term
tmp<fvScalarMatrix> vegetationModel::Sh(volVectorField& U, volScalarField& T) const
{
    Info << "Qs : " << max(Qs(U,T)) << endl;
    return fvm::SuSp((Qs(U,T)/(rhoa_*cpa_))/T,T);
}

// return momentum source term
tmp<fvVectorMatrix> vegetationModel::Su(volVectorField& U) const
{
    return fvm::SuSp(-0.5*Cdf_*LAD_*mag(U), U);
}

// return specific humidity source term
tmp<fvScalarMatrix> vegetationModel::Sq(volVectorField& U, volScalarField& T, volScalarField& q) const
{
    //return E(U,T,q);
    return fvm::SuSp((E(U,T,q)/rhoa_)/q,q);
}

// -----------------------------------------------------------------------------

void vegetationModel::update_Tl(volVectorField& U, volScalarField& T, volScalarField& q)
{
    // net radiation
    volScalarField tmp_Rn("tmp_Rn", Rn(U));

    // vegetation resistance
    volScalarField tmp_ra("tmp_ra", ra(U));

    // bounding LAD
    volScalarField tmp_LAD("tmp_LAD", LAD_);
    dimensionedScalar LADMin("LADMin", dimensionSet(0,-1,0,0,0,0,0), SMALL);
    bound(tmp_LAD, LADMin);

    // Initial leaf temperature
    forAll(tmp_LAD, cellI)
    {
        if (tmp_LAD[cellI] > 100*SMALL)
            Tl_.internalField()[cellI] = T.internalField()[cellI];
    }

    // Info << "WIP <<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    // solve for leaf temperature
    int i = 0; int maxIter = 10;
    while (i < maxIter)
    {
        // Calculate latent heat flux.
        volScalarField tmp_Ql("tmp_Ql", Ql(U, T, q));

        //Calculate new leaf temperature
        volScalarField new_Tl("new_Tl", T + (tmp_Rn - tmp_Ql)*(tmp_ra/(2.0*rhoa_*cpa_*tmp_LAD)) );

        // Determine errors
        scalar maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
        scalar maxRelError = maxError/gMax(mag(new_Tl.internalField()));
        // Iteration info
        Info << "Calculating leaf temp. Iteration i : " << i++
             << ", max. error: " << maxError
             << ", max. rel. error: " << maxRelError << endl;;

        // update leaf temp.
        //Tl_.internalField() = new_Tl.internalField();
        forAll(tmp_LAD, cellI)
        {
            if (tmp_LAD[cellI] > 100*SMALL)
                Tl_.internalField()[cellI] = new_Tl.internalField()[cellI];
        }

        // convergence check
        if (maxRelError < 1e-8)
        {
            Info << "Leaf temperature calculation CONVERGED !!"
                 << ", max. rel. error: " << maxRelError << endl;
            break;
        }

    }

    if (i == maxIter)
    {
        Info << "WARNING!! Leaf temperature not converged !!" << endl;
    }

}


void vegetationModel::testing(volVectorField& U, volScalarField& T, volScalarField& q)
{
    update_Tl(U, T, q);
}

bool vegetationModel::read()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam

// ************************************************************************* //

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
