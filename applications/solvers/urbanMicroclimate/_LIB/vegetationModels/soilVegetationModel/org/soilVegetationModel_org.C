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

#include "soilVegetationModel.H"
#include <ctime>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(soilVegetationModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

soilVegetationModel::soilVegetationModel
(
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& w
	//volScalarField& Tl_
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
    Cf_
    (
        IOobject
        (
            "Cf",
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),
    // pv_
    // (
    //     IOobject
    //     (
    //         "pv",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    // ),
    // pvsat_
    // (
    //     IOobject
    //     (
    //         "pvsat",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    // ),
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
    LAI_
    (
        IOobject
        (
            "LAI",
            "0",//runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    // wsat_
    // (
    //     IOobject
    //     (
    //         "wsat",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0)
    // ),
    Qlat_
    (
        IOobject
        (
            "Qlat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
    // rhosat_
    // (
    //     IOobject
    //     (
    //         "rhosat",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    // ),
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
    Tl_
    (
        IOobject
        (
            "Tl",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,0,1,0,0,0), 0.0)
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
        dimensionedVector("0", dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    )//,
    // VPD_
    // (
    //     IOobject
    //     (
    //         "VPD",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    // )
    {
        // Bounding parameters
		//bound(Tl_, TlMin_);

        Info << " Defined custom vegetation model: soil-vegetation foam" << endl;
    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double soilVegetationModel::calc_pvsat(double& T)
{
    // saturated vapor pressure pws - ASHRAE 1.2
    return exp( - 5.8002206e3/T
                + 1.3914993
                - 4.8640239e-2*T
                + 4.1764768e-5*pow(T,2)
                - 1.4452093e-8*pow(T,3)
                + 6.5459673*log(T) );
}

double soilVegetationModel::calc_pv(double& p, double& w)
{
    // vapor pressure pv
    return (p*w)/(0.6219438209 + w);
}


// calc saturated density of water vapour
double soilVegetationModel::calc_rhosat(double& T)
{
    return calc_pvsat(T)/(461.5*T);
}


// solve radiation
void soilVegetationModel::radiation()
{
    //const fvMesh& vegiMesh =
	//    	mesh_.time().lookupObject<fvMesh>("vegetation");

    //const label patchi = vegiMesh.boundaryMesh().findPatchID("air_to_vegetation");

    //const fvPatch& vegiPatch = vegiMesh.boundary()[patchi];

    //scalarField vegiPatchQr = vegiPatch.lookupPatchField<volScalarField, scalar>("Qr");
    //scalarField vegiPatchQs = vegiPatch.lookupPatchField<volScalarField, scalar>("Qs");
    //scalar integrateQr = gSum(vegiPatch.magSf() * vegiPatchQr);
    //scalar integrateQs = gSum(vegiPatch.magSf() * vegiPatchQs);

	//Info << "test: integrateQr: " << integrateQr << endl;
 	//Info << "test: integrateQs: " << integrateQs << endl;
    //scalar vegiVolume = gSum(pos(Cf_.internalField() - 10*SMALL)*mesh_.V());
	//Info << "test: vegiVolume: " << vegiVolume << endl;

    label timestepsInADay_ = divqrsw.size(); //readLabel(coeffs_.lookup("timestepsInADay"));

    Time& time = const_cast<Time&>(mesh_.time());
    Info << "time.value(): " << time.value();

    label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
    Info << ", 1 timestep: " << timestep;
    timestep = timestep % timestepsInADay_;
    Info << ", 2 timestep: " << timestep << endl;

    //vector sunPos = sunPosVector[timestep];

    scalarList divqrswi =  divqrsw[timestep];

    // radiation density inside vegetation
    forAll(Cf_, cellI)
        if (Cf_[cellI] > 10*SMALL)
          //Rn_[cellI] = -divqrswi[cellI] + (integrateQr)/(vegiVolume);
          Rn_[cellI] = -divqrswi[cellI];// + (integrateQr)/(vegiVolume);

    //Rn_[cellI] = (integrateQr + integrateQs)/(vegiVolume);
    Rn_.correctBoundaryConditions();
    Rn_.write();

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
void soilVegetationModel::resistance(volScalarField& magU, volScalarField& T, volScalarField& w, volScalarField& Tl_)
{
    const double p = 101325;

    double pv, pvsat, VPD;

    // Calculate magnitude of velocity and bounding above Umin
    forAll(Cf_, cellI)
    {
        if (Cf_[cellI] > 10*SMALL)
        {
            //Aerodynamic resistance
            // ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);
            ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

            // Calculate vapor pressure of air
            //ev_[cellI] = q[cellI]*rhoa_.value()*T[cellI]*461.5;
            //pv_[cellI] = p_*w[cellI]/(0.621945+w[cellI]); // orig
            pv = p * w[cellI] / (0.6219438209 + w[cellI]); // new

            // Calculate sat. vapor pressure of air
            //evsat_[cellI] = calc_evsat(T[cellI]); // TODO bug
            //pvsat_[cellI] = calc_pvsat(T[cellI]); // orig
            pvsat = calc_pvsat(T[cellI]); // new

            // Vapor pressure deficit - kPa
            // VPD_[cellI] = (calc_evsat(T[cellI]) - (q[cellI]*rhoa_.value()*T[cellI]*461.5))/1000.0; // kPa
            //VPD_[cellI] = ev_[cellI] - evsat_[cellI];
            //VPD_[cellI] = pvsat_[cellI] - pv_[cellI]; // orig
            VPD = pvsat - pv; // new


            // Stomatal resistance - type 1
            // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
            //rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);
            // rs_[cellI] = rsMin_.value()*(31.0 + Rg_[cellI].component(2))*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rg_[cellI].component(2));

            /*
            // Stomatal resistance - type 2
            // rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2)); // type 2
            //if ((VPD_[cellI]/1000.0) < D0_.value())
            //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()));
            //else
            //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));
            if ((VPD_[cellI]/1000.0) < D0_.value())
                rs_[cellI] = rsMin_.value();//%*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])));
            else
                rs_[cellI] = rsMin_.value();//%*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));

            */
            
            // Stomatal resistance
            rs_[cellI] = rsMin_.value(); // (R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));
            //rs_[cellI] = rsMin*(R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));

        }
    }
    // pv_.correctBoundaryConditions();
    // pvsat_.correctBoundaryConditions();
    // VPD_.correctBoundaryConditions();
    ra_.correctBoundaryConditions();
    rs_.correctBoundaryConditions();
}

// solve vegetation model
//void soilVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& w, volScalarField& Tl_)
void soilVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& w)
{
    // solve radiation within vegetation
    radiation();

    const double p = 101325;
    double pvsat, wsat; // VPD, rhosat, 

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
        // Solve aerodynamic, stomatal resistance
        //std::clock_t startTime = std::clock();
        resistance(magU, T, w, new_Tl);
        //Info << "It took " << (std::clock() - startTime ) / (double) CLOCKS_PER_SEC << " second(s)." << endl;

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                if (i==1)
                    Tl_[cellI] = T[cellI];//*0. + 300.;//T[cellI];

                // Calculate saturated density, specific humidity
                // rhosat_[cellI] = calc_rhosat(Tl_[cellI]); // orig
                //rhosat = calc_rhosat(Tl_[cellI]); // new

                //qsat_[cellI]   = rhosat_[cellI]/rhoa_.value();
                // pvsat_[cellI] = calc_pvsat(Tl_[cellI]); // orig
                pvsat = calc_pvsat(Tl_[cellI]); // new

                //wsat_[cellI] = 0.621945*(pvsat_[cellI]/(p_-pvsat_[cellI])); // ASHRAE 1, eq.23 // orig
                wsat = 0.621945*(pvsat / (p -pvsat) ); // ASHRAE 1, eq.23 // new

                // Calculate transpiration rate
                // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());
                // E_[cellI] = 2.0*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);

				// if (runTime_.value() >= 16800 && runTime_.value() <= 72600) //timesteps are hardcoded - ayk
				// {
                	//E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(wsat_[cellI]-w[cellI])/(ra_[cellI]+rs_[cellI]);// orig
                    E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(wsat-w[cellI])/(ra_[cellI]+rs_[cellI]); // new
                // }
                // else
                // {
                // 	E_[cellI] = 0.0; // No evapotranspiration
                // }
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
            //rhosat_[cellI] = calc_rhosat(Tl_[cellI]); // orig
            //rhosat = calc_rhosat(Tl_[cellI]); // new

            // qsat_[cellI] = rhosat_[cellI]/rhoa_.value();
            // pvsat_[cellI] = calc_pvsat(Tl_[cellI]); // orig
            pvsat = calc_pvsat(Tl_[cellI]); // new
            
            //wsat_[cellI] = 0.621945*(pvsat_[cellI]/(p_-pvsat_[cellI])); // ASHRAE 1, eq.23 // orig
            wsat = 0.621945*(pvsat / (p-pvsat) ); // ASHRAE 1, eq.23 // new

            // Calculate transpiration rate
            // E_[cellI] = LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_.value());
            // E_[cellI] = 2.0*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);

				// if (runTime_.value() >= 16800 && runTime_.value() <= 72600) //timesteps are hardcoded - ayk
				// {
                	// E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(wsat_[cellI]-w[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side
                    E_[cellI] = nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(wsat-w[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side
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
    //rhosat_.correctBoundaryConditions();
    //wsat_.correctBoundaryConditions();
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
tmp<volScalarField> soilVegetationModel::Sh()
{
    Sh_ = Qsen_;
    Sh_.correctBoundaryConditions();
    return Sh_;
}


// return specific humidity source term
tmp<volScalarField> soilVegetationModel::Sw()
{
    Sw_ = E_;
    Sw_.correctBoundaryConditions();
    return Sw_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> soilVegetationModel::Su(volScalarField& rho, volVectorField& U)
{
    return fvm::SuSp(-Cf_*rho*mag(U), U);
}

// tmp<fvScalarMatrix> soilVegetationModel::Sws(volScalarField& rho, volVectorField& U)
// {
//     return E_;
// }

// tmp<volScalarField> soilVegetationModel::Sws()
// {
//     Seta_ = (mtrans()/1000.0)*(RAD*D_eta)/fvc::domainIntegrate(RAD*D_eta);
//     Seta_.correctBoundaryConditions();
//     return Seta_;
// }


// -----------------------------------------------------------------------------

bool soilVegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam



// dimensionedScalar soilVegetationModel::mtrans()
// {
    
//     // Info<< "    Integral of " << E_.name()
//     //     << " over full fluid volume = "
//     //     << gSum(mesh_.V()*E_.internalField()) << " "
//     //     << E_.dimensions()*dimVolume
//     //     << nl;
    
//     // //return gSum(mesh_.V()*E_.internalField());
//     // return dimensionedScalar
//     // (
//     //     "mtrans",
//     //     E_.dimensions()*dimVolume,
//     //     gSum(mesh_.V()*E_.internalField())
//     // );
//     return fvc::domainIntegrate(E_);

// }