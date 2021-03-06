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
    const volScalarField& q,
	        volScalarField& Tl_
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
    beta_
    (
        vegetationProperties_.lookupOrDefault("beta", 0.5)
    ),
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
    l_
    (
        vegetationProperties_.lookup("l")
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
    divqrsw
    (
        IOobject
        (
            "divqrsw",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    cd_
    (
        IOobject
        (
            "cd",
            "0",//runTime_.timeName(),
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),
    LAD_
    (
        IOobject
        (
            "LAD",
            "0",//runTime_.timeName(),
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    Qlat_
    (
        IOobject
        (
            "Qlat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    Qsen_
    (
        IOobject
        (
            "Qsen",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,0,0,0,0,0), 0.0)
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
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
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
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
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
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    ),
    VPD_
    (
        IOobject
        (
            "VPD",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    )
    {
        // Bounding parameters
		    //bound(Tl_, TlMin_);

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
    const fvMesh& vegetationMesh =
	    	mesh_.time().lookupObject<fvMesh>("vegetation");

    const label patchi = vegetationMesh.boundaryMesh().findPatchID("air_to_vegetation");

    const fvPatch& vegiPatch = vegetationMesh.boundary()[patchi];

    scalarField vegiPatchQr = vegiPatch.lookupPatchField<volScalarField, scalar>("Qr");
    scalarField vegiPatchQs = vegiPatch.lookupPatchField<volScalarField, scalar>("Qs");
    scalar integrateQr = gSum(vegiPatch.magSf() * vegiPatchQr);
    scalar integrateQs = gSum(vegiPatch.magSf() * vegiPatchQs);

	   Info << "test: integrateQr: " << integrateQr << endl;
 	   Info << "test: integrateQs: " << integrateQs << endl;
     scalar vegiVolume = gSum(pos(LAD_.internalField() - 10*SMALL)*mesh_.V());
	   Info << "test: vegiVolume: " << vegiVolume << endl;

     label timestepsInADay_ = divqrsw.size(); //readLabel(coeffs_.lookup("timestepsInADay"));

     Time& time = const_cast<Time&>(mesh_.time());
     Info << "time.value(): " << time.value();

     label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
     Info << ", 1 timestep: " << timestep;
     timestep = timestep % timestepsInADay_;
     Info << ", 2 timestep: " << timestep << endl;

     //vector sunPos = sunPosVector[timestep];

     scalarList divqrswi =  divqrsw[timestep];

     //Info << "I0 " << -sunPos*IDN[timestep] << endl;

    //  forAll(LAD_, cellI)
    //  {
    //     LAI_[cellI] = LAIList[timestep][cellI];
    //     Rg_[cellI] = -sunPos*IDN[timestep]*Foam::exp(-beta_*LAI_[cellI]);
    //  }
    //
    //  Rg_.correctBoundaryConditions();
    //  Rg_.write();
    //
    //  volScalarField divRg(fvc::div(Rg_));
    //  divRg.correctBoundaryConditions();
    //  divRg.write();




     //Rg_.correctBoundaryConditions();

     //volTensorField gradRg = fvc::grad(Rg_);
     //gradRg.correctBoundaryConditions();

     // radiation density inside vegetation
     forAll(LAD_, cellI)
     {
        if (LAD_[cellI] > 10*SMALL)
        {
            //Rn_[cellI] = - divqrswi[cellI] + (integrateQr + integrateQs)/(vegiVolume); // direct and diffuse
            Rn_[cellI] = (integrateQr + integrateQs)/(vegiVolume); // direct and diffuse
            //Rn_[cellI] = - divqrswi[cellI] + (integrateQs)/(vegiVolume); // direct and diffuse

        }

     }

     Rn_.correctBoundaryConditions();

     Info << "gMax(Rn) = " << gMax(Rn_) << ", gMin(Rn) = " << gMin(Rn_) << endl;
     Info << "sum qrswi = " << gSum(divqrswi*mesh_.V()) << endl;

     //Rn_.write();
     //LAI_.write();

	//Info << vegiPatch.Cf() << endl;
	// Now I have the face centres for air_to_vegetation boundary.

	/*forAll(Cf_, cellI)
	{
        if (Cf_[cellI] > 10*SMALL)
        {
        	Info << "Cf_[cellI]: " << Cf_[cellI];
        	forAll(vegiPatch, faceI)
        	{
		        if (Cf_[cellI] > 10*SMALL)
        		{
					check normal vectors at each vegipatch face and see if they are pos or neg to understand if top boundary or bottom, or left or right
        		}
        	}
        }
	}*/

}

// solve aerodynamic resistance
void vegetationModel::resistance(volScalarField& magU, volScalarField& T, volScalarField& q, volScalarField& Tl_)
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
void vegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& q, volScalarField& Tl_)
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

				        //if (runTime_.value() >= 16800 && runTime_.value() <= 72600) //timesteps are hardcoded - ayk
				        //{
                	E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);
                //}
                //else
                //{
                //	E_[cellI] = 0.0; // No evapotranspiration
                //}
                //E_[cellI] = 0.0; // No evapotranspiration

                // Calculate latent heat flux
                //Qlat_[cellI] = lambda_.value()*E_[cellI];
                Qlat_[cellI] = lambda_.value()*E_[cellI];

                // Calculate new leaf temperature
                //new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Qlat_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));
                new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Qlat_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));

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
    //Tl_.write();

    // Iteration info
    Info << "Vegetation model:  Solving for Tl, Final residual = " << maxError
         << ", Final relative residual = " << maxRelError
         << ", No Iterations " << i << endl;

    Info << "temperature parameters: min Tl = " << gMin(Tl_) << ", max Tl = " << gMax(Tl_)
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

				// if (runTime_.value() >= 16800 && runTime_.value() <= 72600) //timesteps are hardcoded - ayk
				// {
      	     E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side
        // }
        // else
        // {
        // 	E_[cellI] = 0.0; // No evapotranspiration
        // }

            //E_[cellI] = 0.0; // no evapotranspiration
            // TODO: flag for no transpiration, one side, both side

            // Calculate latent heat flux
            Qlat_[cellI] = lambda_.value()*E_[cellI];

            // Calculate sensible heat flux
            Qsen_[cellI] = 2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]*(Tl_[cellI]-T[cellI])/ra_[cellI];
        }
    }
    rhosat_.correctBoundaryConditions();
    qsat_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    Qlat_.correctBoundaryConditions();
    Qsen_.correctBoundaryConditions();

    // Iteration info
    // Info << "              Vegetation model:  max. Rn = " << max(mag(Rn_))
    //      << "; max. Qlat = " << max(mag(Qlat_))
    //      << "; max. Qsen = " << max(mag(Qsen_))
    //      << "; error: max. Esum = " << max(mag(Rn_.internalField() - Qsen_.internalField()- Qlat_.internalField())) << endl;

}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> vegetationModel::Sh()
{
    Sh_ = Qsen_;
    Sh_.correctBoundaryConditions();
    return Sh_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> vegetationModel::Su(volScalarField& rho, volVectorField& U)
{
    Su_ = -cd_*LAD_*rho*mag(U)*U;
    Su_.correctBoundaryConditions();
    return fvm::SuSp(-cd_*LAD_*rho*mag(U), U);
}

// return specific humidity source term
tmp<volScalarField> vegetationModel::Sw()
{
    Sw_ = E_;
    Sw_.correctBoundaryConditions();
    return Sw_;
}

// -----------------------------------------------------------------------------

bool vegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam
