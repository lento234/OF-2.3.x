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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(soilVegetationModel, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

double soilVegetationModel::calc_pvsat(const double& T)
{
    // saturated vapor pressure pws - ASHRAE 1.2
    return exp( - 5.8002206e3/T
                + 1.3914993
                - 4.8640239e-2*T
                + 4.1764768e-5*pow(T,2)
                - 1.4452093e-8*pow(T,3)
                + 6.5459673*log(T) );
}

double soilVegetationModel::calc_pv(const double& p, const double& w)
{
    // vapor pressure pv
    return (p*w)/(0.6219438209 + w);
}

double soilVegetationModel::calc_rhosat(const double& T)
{
    return calc_pvsat(T)/(461.5*T);
}

// void soilVegetationModel::calc_radiation()
// {
//     label timestepsInADay_ = divqrsw.size();

//     Time& time = const_cast<Time&>(mesh_.time());
//     Info << "Vegetation: [Radiation] :: time.value(): " << time.value();

//     label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
//     Info << ", 1 timestep: " << timestep;
//     timestep = timestep % timestepsInADay_;
//     Info << ", 2 timestep: " << timestep << endl;

//     scalarList divqrswi =  divqrsw[timestep];

//     // radiation density inside vegetation
//     forAll(Cf_, cellI)
//         if (Cf_[cellI] > 10*SMALL)
//             qrad_leaf_[cellI] = (- divqrswi[cellI] )/ LAD_[cellI];

//     qrad_leaf_.correctBoundaryConditions();

// }


void soilVegetationModel::calc_radiation()
{
    // ---------------------------------------------
    // Direct+Diffused long-wave and short-wave
    // ---------------------------------------------
    
    // Direct+Diffused long-wave and short-wave radiative heat flux at vegetation boundary
    const fvMesh& vegiMesh = mesh_.time().lookupObject<fvMesh>("vegetation"); // "vegetation" a.k.a radiation domain mesh
    const label patchi = vegiMesh.boundaryMesh().findPatchID("air_to_vegetation"); // patch id
    const fvPatch& vegiPatch = vegiMesh.boundary()[patchi]; // patch boundary of vegetation
    scalarField vegiPatchQr = vegiPatch.lookupPatchField<volScalarField, scalar>("Qr"); // long-wave (diffused)
    scalarField vegiPatchQs = vegiPatch.lookupPatchField<volScalarField, scalar>("Qs"); // short-wave (diffused + direct)

    // Net radiation inside vegetation : surface integral
    scalar integrateQr = gSum(vegiPatch.magSf() * vegiPatchQr);
    scalar integrateQs = gSum(vegiPatch.magSf() * vegiPatchQs);
    scalar vegiVolume  = gSum(pos(Cf_.internalField() - 10*SMALL)*mesh_.V());

    
    Info << "Vegetation: [Radiation] :: int. Qr = " << integrateQr
         << ", int. Qs = " << integrateQs 
         << ", vol. = " << vegiVolume; 

    if ((integrateQs/vegiVolume) > (100 * SMALL))
    {
        isDayTime_ = true;
        Info << ", it is day time.";
    }
    else
    {
        isDayTime_ = false;
        Info << ", it is night time.";
    }
    Info << endl;
    
    // ---------------------------------------------
    // Direct short-wave radiation through beers law
    // ---------------------------------------------

    // Time of the day
    label timestepsInADay_ = divqrsw.size();
    Time& time = const_cast<Time&>(mesh_.time());
    Info << "Vegetation: [Radiation] :: time.value() = " << time.value();
    label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
    Info << ", 1 timestep = " << timestep;
    timestep = timestep % timestepsInADay_;
    Info << ", 2 timestep = " << timestep << endl;

    // Net radiation absorbed by vegetation
    scalarList divqrswi =  divqrsw[timestep];

    // radiation density inside vegetation
    forAll(Cf_, cellI)
        if (Cf_[cellI] > 10*SMALL)
            qrad_leaf_[cellI] = (- divqrswi[cellI] )/ LAD_[cellI];
            //qrad_leaf_[cellI] =  ( -divqrswi[cellI]  + (integrateQr/vegiVolume) ) /  LAD_[cellI];
    qrad_leaf_.correctBoundaryConditions();

    // Integrate direct short-wave radiation absorbed by vegetation
    dimensionedScalar integrateQrsw = fvc::domainIntegrate(qrad_leaf_ * LAD_);

    
    

    Info << "Vegetation: [Radiation] :: int. Qrsw = " << integrateQrsw << endl;

    // Missing short-wave radiation, diffused and otherwise.
    scalar missingQs = integrateQs - integrateQrsw.value(); // needs to be added to satisfy energy balance, sadly

    forAll(Cf_, cellI)
        if (Cf_[cellI] > 10*SMALL)
            qrad_leaf_[cellI] += ((integrateQr + missingQs)/vegiVolume)/ LAD_[cellI];
            //qrad_leaf_[cellI] =  ( -divqrswi[cellI]  + (integrateQr/vegiVolume) ) /  LAD_[cellI];
    qrad_leaf_.correctBoundaryConditions();

    Info << "Vegetation: [Radiation] :: update. int. Qrsw = " << fvc::domainIntegrate(qrad_leaf_ * LAD_) 
         << ", Qr + Qs = " << integrateQr + integrateQs << endl;

    // update internal clock
    internalTime = time.value();

}

/*
double soilVegetationModel::calc_resistance_aerodynamic(const double& magU)
{
    return (C_.value() * pow(l_.value() / magU, 0.5));
}
*/

// 
double soilVegetationModel::calc_conductance_aerodynamic(const double& magU)
{
    return max( ( pow(magU / l_.value(), 0.5) / C_.value() ), SMALL);
}


double soilVegetationModel::calc_conductance_stomatal(const double& pv, const double& pvsat, const double& T, const int& cellI)
{
   
    // Vapor pressure deficit
    //VPD = pvsat - pv; // new

    // Stomatal resistance - type 1
    // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
    // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);
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

    if (isDayTime_)
    {
        // constant Stomatal resistance
        // case STOMATAL_RESISTANCE_CONSTANT: return rsMin_.value();
        return 1.0/rsMin_.value();    
    }
    else
    {
        return SMALL;
    }
    
    // Bruse and Fleer (1998)
    // return rsMin_.value()* // (R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));
    
}

// double soilVegetationModel::calc_resistance_stomatal(const double& pv, const double& pvsat, const double& T, const int& cellI)
// {
   
//     // Vapor pressure deficit
//     //VPD = pvsat - pv; // new

//     // Stomatal resistance - type 1
//     // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
//     // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);
//     // rs_[cellI] = rsMin_.value()*(31.0 + Rg_[cellI].component(2))*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rg_[cellI].component(2));

//     /*
//     // Stomatal resistance - type 2
//     // rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2)); // type 2
//     //if ((VPD_[cellI]/1000.0) < D0_.value())
//     //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()));
//     //else
//     //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));
//     if ((VPD_[cellI]/1000.0) < D0_.value())
//         rs_[cellI] = rsMin_.value();//%*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])));
//     else
//         rs_[cellI] = rsMin_.value();//%*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));

//     */

    
//     // constant Stomatal resistance
//     // case STOMATAL_RESISTANCE_CONSTANT: return rsMin_.value();
//     return rsMin_.value();    
    
//     // Bruse and Fleer (1998)
//     // return rsMin_.value()* // (R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));
    
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

soilVegetationModel::soilVegetationModel
(
    const volVectorField& U,
    const rhoThermo& thermo,
    const volScalarField& w,
    const volScalarField& Ts, // soil added
    const volScalarField& ws // soil added
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
    runTimeSoil_(Ts.time()),
    mesh_(U.mesh()),
    meshSoil_(Ts.mesh()),
    internalTime(-1.0),
    isDayTime_(false),
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
    //TlMin_("TlMin", dimTemperature, SMALL),
    UMin_("UMin", dimVelocity, SMALL),
    lambda_
    (
        vegetationProperties_.lookup("lambda")
    ),
    mtrans_
    (
        "mtrans",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
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
            runTime_.timeName(),//"0",//runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gv_leaf_
    (
        IOobject
        (
            "gv_leaf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
    ),
    LAD_
    (
        IOobject
        (
            "LAD",
            runTime_.timeName(),//"0",//runTime_.timeName(),
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
            runTime_.timeName(),//"0",//runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qlat_leaf_
    (
        IOobject
        (
            "qlat_leaf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
    ),
    qsen_leaf_
    (
        IOobject
        (
            "qsen_leaf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
    ),
    ga_
    (
        IOobject
        (
            "ga",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    gs_
    (
        IOobject
        (
            "gs",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    qrad_leaf_
    (
        IOobject
        (
            "qrad_leaf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
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
    h_ch_
    (
        IOobject
        (
            "h_ch",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,0,-3,-1,0,0,0), 0.0)
    ),
    h_cm_
    (
        IOobject
        (
            "h_cm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    Sh_
    (
        IOobject
        (
            "Sh",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE//AUTO_WRITE
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
            IOobject::NO_WRITE//AUTO_WRITE
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
            IOobject::NO_WRITE//AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    ),
    Sws_
    (
        IOobject
        (
            "Sws",
            runTimeSoil_.timeName(),
            meshSoil_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshSoil_,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),
    RAD_
    (
        IOobject
        (
            "RAD",
            runTimeSoil_.timeName(),
            meshSoil_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshSoil_
    )     
    {
        Info << " Defined custom vegetation model: soil-vegetation foam" << endl;
    }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// solve vegetation model
//void soilVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& w, volScalarField& Tl_)
//void soilVegetationModel::solve(volVectorField& U, volScalarField& T, volScalarField& w)
void soilVegetationModel::solve(const volVectorField& U, const rhoThermo& thermo, const volScalarField& w)
{
    const volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& rho = thermo.rho();
    double pvsat, pv, gnet;
    double Ra = 287.042;
    double Rv = 461.524;

    // Initialize radiation and leaf temperature if internalTime is changed
    if (mesh_.time().value() != internalTime)
    {
        Info << "Vegetation: [Radiation] :: updating radiation field" 
             << ", internal time = " << internalTime << endl;
        calc_radiation();

        // Initialize leaf temperature
        forAll(Tl_, cellI)
            if (LAD_[cellI] > 10*SMALL)
                Tl_[cellI] = T[cellI];
    }    

    Info << "Vegetation: [Leaf Energy Balance] :: updating leaf temperature "
         << ", internal time = " << internalTime << endl;

    // Magnitude of velocity
    volScalarField magU("magU", mag(U));

    // Bounding velocity
    // bound(magU, UMin_); // no longer needed

    // solve aerodynamic, stomatal resistance
    //resistance(U,T);
    volScalarField new_Tl("new_Tl", Tl_);

    scalar maxError, maxRelError;
    int i;

    // solve leaf temperature, iteratively.
    int maxIter = 500;
    for (i=1; i<=maxIter; i++)
    {
        // Solve aerodynamic, stomatal resistance
        //std::clock_t startTime = std::clock();
        //resistance(magU, T, w, new_Tl);
        //Info << "It took " << (std::clock() - startTime ) / (double) CLOCKS_PER_SEC << " second(s)." << endl;

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                // if (i==1)
                //     Tl_[cellI] = T[cellI];//*0. + 300.;//T[cellI];

                // vapour pressure
                pv = p[cellI] * w[cellI] / (Ra/Rv + w[cellI]); // vapour pressure

                // saturation vapour pressure
                pvsat = calc_pvsat(Tl_[cellI]);

                // aerodynamic resistance
                //ra_[cellI] = calc_resistance_aerodynamic(magU[cellI]);
                // Aerodynamic conductance (boundary-layer conductance) m/s
                ga_[cellI] = calc_conductance_aerodynamic(magU[cellI]);

                // stomatal resistance
                // rs_[cellI] = calc_resistance_stomatal(pv, pvsat, T[cellI], cellI);
                // stomatal conductance m/s
                gs_[cellI] = calc_conductance_stomatal(pv,pvsat,T[cellI],cellI);

                // net heat/vapor conductance m/s
                gnet = (gs_[cellI] * ga_[cellI]) / (gs_[cellI] + ga_[cellI]) ;

                // convective heat transfer coefficient
                //h_ch_[cellI] = (2.0*rho[cellI]*cpa_.value())/ra_[cellI];
                h_ch_[cellI] = 2.0 * ga_[cellI] * rho[cellI] * cpa_.value();

                // convective mass transfer coefficient
                //h_cm_[cellI] = (rho[cellI]*Ra)/(p[cellI]*Rv*(ra_[cellI]+rs_[cellI]));
                h_cm_[cellI] = gnet * (rho[cellI]*Ra / (p[cellI]*Rv) ); //*(ra_[cellI]+rs_[cellI]));

                // Calculate transpiration rate (mass flux rate)
                gv_leaf_[cellI] = nEvapSides_.value() * h_cm_[cellI] * (pvsat - pv);
                
                // Calculate latent heat flux
                qlat_leaf_[cellI] = lambda_.value() * gv_leaf_[cellI];

                // Calculate sensible heat flux
                qsen_leaf_[cellI] = h_ch_[cellI] * (Tl_[cellI] - T[cellI]);

                // Calculate new leaf temperature
                //new_Tl[cellI] = T[cellI] + (qrad_leaf_[cellI] - qlat_leaf_[cellI]) * (ra_[cellI]/(2.0*rho[cellI]*cpa_.value()));
                new_Tl[cellI] = T[cellI] + (qrad_leaf_[cellI] - qlat_leaf_[cellI]) / h_ch_[cellI];

            }
        }
        
        // Check rel. L-infinity error
        maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
        maxRelError = maxError/gMax(mag(new_Tl.internalField()));

        // Info
        Info << "Vegetation: [Leaf Energy Balance] :: Iteration = " << i
             << " max. Tl = " << gMax(new_Tl)
             << ", error Tl = "   << maxError
             << " LEB error = " << gMax(qrad_leaf_.internalField() - qsen_leaf_.internalField() - qlat_leaf_.internalField())
             << endl;
        
        // convergence check
        if ((maxRelError < 1e-8) && (maxError < 1e-8))
            break;
        else
            forAll(Tl_, cellI)
                Tl_[cellI] = 0.5*Tl_[cellI]+0.5*new_Tl[cellI]; // stabilized // update leaf temp.       
    }

    // Correct boundary conditions
    ga_.correctBoundaryConditions();
    gs_.correctBoundaryConditions();
    gv_leaf_.correctBoundaryConditions();
    qlat_leaf_.correctBoundaryConditions();
    qsen_leaf_.correctBoundaryConditions();

    // Iteration info
    Info << "Vegetation: [Leaf Energy Balance] :: Final residual = " << maxError
         << ", Final rel. residual = " << maxRelError
         << ", No Iterations " << i 
         << endl;
    
    Info << "Vegetation: [Leaf Energy Balance] :: int. a*qr = " << fvc::domainIntegrate(qrad_leaf_ * LAD_).value()
         << ", int. a*qs = " << fvc::domainIntegrate(qsen_leaf_ * LAD_).value()
         << ", int. a*ql = " << fvc::domainIntegrate(qlat_leaf_ * LAD_).value() 
         << endl;

}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> soilVegetationModel::Sh()
{
    //Sh_ = Qsen_;
    Sh_ = LAD_ * qsen_leaf_;
    Sh_.correctBoundaryConditions();
    return Sh_;
}


// return humidity source term (moisture content)
tmp<volScalarField> soilVegetationModel::Sw()
{
    //Sw_ = E_;
    Sw_ = LAD_ * gv_leaf_;
    Sw_.correctBoundaryConditions();

    return Sw_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> soilVegetationModel::Su(volScalarField& rho, volVectorField& U)
{
    // Su_ = -Cf_*rho*mag(U)*U;
    // Su_.correctBoundaryConditions();
    return fvm::SuSp(-Cf_*rho*mag(U), U); // Cf_ = cd*LAD_
}


tmp<volScalarField> soilVegetationModel::Sws(volScalarField& Kl, volScalarField& Cl)
{
    // Hydraulic liquid diffusivity, m2/s
    volScalarField Dl("Dl", Kl/Cl);

    // Net transpiration rate kg/s
    dimensionedScalar mtrans("mtrans", fvc::domainIntegrate(LAD_ * gv_leaf_));
 
    // source term for water uptake due to roots kg/(m3s)
    Sws_ = - (mtrans * RAD_ * Dl) / fvc::domainIntegrate(RAD_ * Dl);
    
    Info << "Vegetation: [soil root uptake] :: Sws_ max = " << gMax(Sws_) 
         << ", min: " << gMin(Sws_) << endl;

    return Sws_;
}


// -----------------------------------------------------------------------------

bool soilVegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam




// void soilVegetationModel::calc_radiation()
// {
//     //const fvMesh& vegiMesh =
// 	//    	mesh_.time().lookupObject<fvMesh>("vegetation");

//     //const label patchi = vegiMesh.boundaryMesh().findPatchID("air_to_vegetation");

//     //const fvPatch& vegiPatch = vegiMesh.boundary()[patchi];

//     //scalarField vegiPatchQr = vegiPatch.lookupPatchField<volScalarField, scalar>("Qr");
//     //scalarField vegiPatchQs = vegiPatch.lookupPatchField<volScalarField, scalar>("Qs");
//     //scalar integrateQr = gSum(vegiPatch.magSf() * vegiPatchQr);
//     //scalar integrateQs = gSum(vegiPatch.magSf() * vegiPatchQs);

// 	//Info << "test: integrateQr: " << integrateQr << endl;
//  	//Info << "test: integrateQs: " << integrateQs << endl;
//     //scalar vegiVolume = gSum(pos(Cf_.internalField() - 10*SMALL)*mesh_.V());
// 	//Info << "test: vegiVolume: " << vegiVolume << endl;

//     label timestepsInADay_ = divqrsw.size(); //readLabel(coeffs_.lookup("timestepsInADay"));

//     Time& time = const_cast<Time&>(mesh_.time());
//     Info << "time.value(): " << time.value();

//     label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
//     Info << ", 1 timestep: " << timestep;
//     timestep = timestep % timestepsInADay_;
//     Info << ", 2 timestep: " << timestep << endl;

//     //vector sunPos = sunPosVector[timestep];

//     scalarList divqrswi =  divqrsw[timestep];

//     // radiation density inside vegetation
//     forAll(Cf_, cellI)
//         if (Cf_[cellI] > 10*SMALL)
//           //Rn_[cellI] = -divqrswi[cellI] + (integrateQr)/(vegiVolume);
//           //Rn_[cellI] = -divqrswi[cellI] / LAD_;// + (integrateQr)/(vegiVolume);
//           qrad_leaf_[cellI] = (- divqrswi[cellI] )/ LAD_[cellI];

//     //Rn_[cellI] = (integrateQr + integrateQs)/(vegiVolume);
//     // Rn_.correctBoundaryConditions();
//     // Rn_.write();
//     qrad_leaf_.correctBoundaryConditions();

// 	//Info << vegiPatch.Cf() << endl;
// 	// Now I have the face centres for air_to_vegetation boundary.

// 	/*forAll(Cf_, cellI)
// 	{
//         if (Cf_[cellI] > 10*SMALL)
//         {
//         	Info << "Cf_[cellI]: " << Cf_[cellI];
//         	forAll(vegiPatch, faceI)
//         	{
// 		        if (Cf_[cellI] > 10*SMALL)
//         		{
// 					check normal vectors at each vegipatch face and see if they are pos or neg to understand if top boundary or bottom, or left or right
//         		}
//         	}
//         }
// 	}*/

// }



// // solve aerodynamic resistance
// void soilVegetationModel::resistance(const volScalarField& magU, const volScalarField& p, const volScalarField& T, const volScalarField& w, const volScalarField& Tl_)
// {
//     //const double p = 101325;

//     double pv, pvsat, VPD;
//     double Ra = 287.042;
//     double Rv = 461.524;

//     // Calculate magnitude of velocity and bounding above Umin
//     forAll(Cf_, cellI)
//     {
//         if (Cf_[cellI] > 10*SMALL)
//         {
//             //Aerodynamic resistance
//             // ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);
//             ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

//             // Calculate vapor pressure of air
//             //ev_[cellI] = q[cellI]*rhoa_.value()*T[cellI]*461.5;
//             //pv_[cellI] = p_*w[cellI]/(0.621945+w[cellI]); // orig
//             pv = p[cellI] * w[cellI] / (Ra/Rv + w[cellI]); // new

//             // Calculate sat. vapor pressure of air
//             //evsat_[cellI] = calc_evsat(T[cellI]); // TODO bug
//             //pvsat_[cellI] = calc_pvsat(T[cellI]); // orig
//             pvsat = calc_pvsat(T[cellI]); // new

//             // Vapor pressure deficit - kPa
//             // VPD_[cellI] = (calc_evsat(T[cellI]) - (q[cellI]*rhoa_.value()*T[cellI]*461.5))/1000.0; // kPa
//             //VPD_[cellI] = ev_[cellI] - evsat_[cellI];
//             //VPD_[cellI] = pvsat_[cellI] - pv_[cellI]; // orig
//             VPD = pvsat - pv; // new


//             // Stomatal resistance - type 1
//             // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
//             //rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);
//             // rs_[cellI] = rsMin_.value()*(31.0 + Rg_[cellI].component(2))*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rg_[cellI].component(2));

//             /*
//             // Stomatal resistance - type 2
//             // rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2)); // type 2
//             //if ((VPD_[cellI]/1000.0) < D0_.value())
//             //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()));
//             //else
//             //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));
//             if ((VPD_[cellI]/1000.0) < D0_.value())
//                 rs_[cellI] = rsMin_.value();//%*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])));
//             else
//                 rs_[cellI] = rsMin_.value();//%*((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));

//             */
            
//             // Stomatal resistance
//             rs_[cellI] = rsMin_.value(); // (R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));
//             //rs_[cellI] = rsMin*(R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));

//         }
//     }
//     // pv_.correctBoundaryConditions();
//     // pvsat_.correctBoundaryConditions();
//     // VPD_.correctBoundaryConditions();
//     ra_.correctBoundaryConditions();
//     rs_.correctBoundaryConditions();
// }

