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


using namespace Foam::constant;
//using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(soilVegetationModel, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class tmpClass>
void soilVegetationModel::assertUnit(const tmpClass& sourceVar, const dimensionSet& targetUnit)
{
    if (sourceVar.dimensions() != targetUnit)
    {
        FatalError
            << "Vegetation : [Status]      :: "
            << sourceVar.name() << " units should be "
            << targetUnit << " and not "
            << sourceVar.dimensions()
            << exit(FatalError);
    }
    else //debug
    {
        Info << sourceVar << endl;
    }
}

void soilVegetationModel::writeVegetationProperties()
{
    // Update dictionary
    //varyingVegetationProperties_.set("c", c_);
    //varyingVegetationProperties_.set("cao", cao_);
    varyingVegetationProperties_.set("psi_L", psi_L_);
    varyingVegetationProperties_.set("psi_R", psi_R_);
    varyingVegetationProperties_.set("lambda_", lambda_);

    //Info
    Info << "Vegetation : [Status]      :: Writing varying vegetation properties" << endl;

    // Export to file
    varyingVegetationProperties_.regIOobject::write();
}

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


double soilVegetationModel::calc_VPD(const double& T, const double& p, const double& w)
{
    // saturation vapour pressure [Pa]
    double pvs = calc_pvsat(T);
    // vapour pressure [Pa]
    double pv = calc_pv(p, w);

    return (pvs - pv)/Pstd_.value();
}

 
double soilVegetationModel::calc_ga(const double& magU, const double& rho)
{
    return (rho/(Mco2_.value()*C_.value())) * pow(magU / l_.value(), 0.5);
}


void soilVegetationModel::calc_soilProperties(const volScalarField& Kl)
{
    //- alpha
    const dimensionedScalar alpha("alpha", pow(Hr_/RAI_, 0.5)/pow(2.0*r_, 0.5));

    // Absolute value of gravity [m/s2]
    const dimensionedScalar gabs("gabs", mag(g_));

    // Calculate conductivies
    double K, ks, meshDL, kr;
    forAll(RAD_, cellI)
    {
        if (RAD_[cellI] > minThreshold)
        {
            //- Hydraulic conductivity [m/s]
            K = Kl[cellI] * gabs.value();

            // Soil-root interface conductance [1/s]
            ks = alpha.value() * K * RAD_[cellI];

            // Characterstic cell size
            meshDL = pow(meshSoil_.V()[cellI],1.0/3.0);

            // Root system condutance (1/s)
            kr = RAD_[cellI] * meshDL / beta_.value();

            // Effective soil-root conductivity (s/m)
            gsr_[cellI] = (ks * kr) / ((ks + kr) * gabs.value());
        }

    }

}


void soilVegetationModel::calc_marginalWUE()
{
    // Average leaf water potential of the last 24 hours
    scalar psi_L24_avg = average(psi_L24_);

    Info << "Vegetation : [MWUE]        :: psi_L24_avg =  " << psi_L24_avg  << endl;

    // Calculate Marginal Water use efficiency (lamba) 
    lambda_.value() = lambda_max_.value() * (c_.value() / c_star_.value()) * exp ( - betaL_.value() * pow(psi_L24_avg - psi_Lmax_.value(), 2));
}


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

    
    Info << "Vegetation : [Radiation]   :: int. Qr = " << integrateQr
         << ", int. Qs = " << integrateQs 
         << ", vol. = " << vegiVolume << endl; 

    // if ((integrateQs/vegiVolume) > (100 * SMALL))
    // {
    //     isDayTime_ = true;
    //     Info << ", it is day time.";
    // }
    // else
    // {
    //     isDayTime_ = false;
    //     Info << ", it is night time.";
    // }
    // Info << endl;
    
    // ---------------------------------------------
    // Direct short-wave radiation through beers law
    // ---------------------------------------------

    // Time of the day
    label timestepsInADay_ = divqrsw.size();
    Time& time = const_cast<Time&>(mesh_.time());
    Info << "Vegetation : [Radiation]   :: time.value() = " << time.value();
    label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
    Info << ", 1 timestep = " << timestep;
    timestep = timestep % timestepsInADay_;
    Info << ", 2 timestep = " << timestep << endl;

    // Net radiation absorbed by vegetation
    scalarList divqrswi =  divqrsw[timestep];

    // radiation density inside vegetation
    forAll(LAD_, cellI)
        if (LAD_[cellI] > minThreshold)
            qrad_leaf_[cellI] = (- divqrswi[cellI] )/ LAD_[cellI];
            //qrad_leaf_[cellI] =  ( -divqrswi[cellI]  + (integrateQr/vegiVolume) ) /  LAD_[cellI];
    qrad_leaf_.correctBoundaryConditions();

    // Integrate direct short-wave radiation absorbed by vegetation
    dimensionedScalar integrateQrsw = fvc::domainIntegrate(qrad_leaf_ * LAD_);


    Info << "Vegetation : [Radiation]   :: int. Qrsw = " << integrateQrsw << endl;

    // Missing short-wave radiation, diffused and otherwise.
    scalar missingQs = integrateQs - integrateQrsw.value(); // needs to be added to satisfy energy balance, sadly

    forAll(LAD_, cellI)
        if (LAD_[cellI] > minThreshold)
            qrad_leaf_[cellI] += ((integrateQr + missingQs)/vegiVolume)/ LAD_[cellI];
            //qrad_leaf_[cellI] =  ( -divqrswi[cellI]  + (integrateQr/vegiVolume) ) /  LAD_[cellI];
    qrad_leaf_.correctBoundaryConditions();

    Info << "Vegetation : [Radiation]   :: update. int. Qrsw = " << fvc::domainIntegrate(qrad_leaf_ * LAD_) 
         << ", Qr + Qs = " << integrateQr + integrateQs << endl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

soilVegetationModel::soilVegetationModel
(
    const volVectorField& U,    //- Fluid velocity [m/s]
    const rhoThermo& thermo,    //- Fluid thermodynamic properties
    const volScalarField& w,    //- Fluid absolute humidity [kg/kg]
    const volScalarField& Ts,   //- Soil temperature [K]
    const volScalarField& ws,   //- Soil moisture [kg/m3]
    const volScalarField& pc    //- Soil capillary pressure [Pa, kg/ms2]
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
    UMin_("UMin", dimVelocity, SMALL),
    minThreshold(10*SMALL),

    rhow_("rhow", dimDensity, 1000.0),
    Pstd_("Pstd", dimPressure, 101300),
    Ra_("Ra", dimGasConstant, 287.042),
    Rv_("Rv", dimGasConstant, 461.524),
    cpa_("cpa", dimSpecificHeatCapacity, 1003.5),
    Mw_("Mw", dimMass/dimMoles, 0.01802),
    Mco2_("Mco2", dimMass/dimMoles, 0.04401),
    Lv_("Lv", dimEnergy/dimMass, 2.5e6),

    c_
    (
        vegetationProperties_.lookup("c")
    ),
    cao_
    (
        vegetationProperties_.lookup("cao")
    ),
    g_
    (
        vegetationProperties_.lookup("g")
    ),
    C_
    (
        vegetationProperties_.lookup("C")
    ),
    l_
    (
        vegetationProperties_.lookup("l")
    ),
    Hr_
    (
        vegetationProperties_.lookup("Hr")
    ),
    RAI_
    (
        vegetationProperties_.lookup("RAI")
    ),
    r_
    (
        vegetationProperties_.lookup("r")
    ),    
    beta_
    (
        vegetationProperties_.lookup("beta")
    ),
    psi_L24_0_
    (
        vegetationProperties_.lookup("psi_L24_0")
    ),
    timestepsInADay_
    (
        vegetationProperties_.lookupOrDefault("timestepsInADay", 24)
    ),
    lambda_max_
    (
        vegetationProperties_.lookup("lambda_max")
    ),    
    c_star_
    (
        vegetationProperties_.lookup("c_star")
    ),    
    betaL_
    (
        vegetationProperties_.lookup("betaL")
    ),  
    psi_Lmax_
    (
        vegetationProperties_.lookup("psi_Lmax")
    ),           
     

    psi_L_("psi_L", dimPressure, 0.0),
    psi_R_("psi_R", dimPressure, 0.0),
    psi_L24_(timestepsInADay_,psi_L24_0_.value()),
    lambda_("lambda", dimMoles/dimMoles, 0.0),

    varyingVegetationProperties_
    (
        IOobject
        (
            "varyingVegetationProperties",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
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

    // isDayTime_(false),
    // a1_
    // (
    //     vegetationProperties_.lookup("a1")
    // ),
    // a2_
    // (
    //     vegetationProperties_.lookup("a2")
    // ),
    // a3_
    // (
    //     vegetationProperties_.lookup("a3")
    // ),
    // cpa_
    // (
    //     vegetationProperties_.lookup("cpa")
    // ),
    
    // D0_
    // (
    //     vegetationProperties_.lookup("D0")
    // ),
    // nEvapSides_
    // (
    //     vegetationProperties_.lookup("nEvapSides")
    // ),
    // H_
    // (
    //     vegetationProperties_.lookup("H")
    // ),
    // kc_
    // (
    //     vegetationProperties_.lookup("kc")
    // ),
    
    // Rg0_
    // (
    //     vegetationProperties_.lookup("Rg0")
    // ),
    // Rl0_
    // (
    //     vegetationProperties_.lookup("Rl0")
    // ),
    // rhoa_
    // (
    //     vegetationProperties_.lookup("rhoa")
    // ),
    // rsMin_
    // (
    //     vegetationProperties_.lookup("rsMin")
    // ),
    
    // lambda_
    // (
    //     vegetationProperties_.lookup("lambda")
    // ),
    // wPWP_
    // (
    //     vegetationProperties_.lookup("wPWP")
    // ),    
    // mtrans_
    // (
    //     "mtrans",
    //     dimensionSet(1,0,-1,0,0,0,0),
    //     0.0
    // ),
    // divqrsw
    // (
    //     IOobject
    //     (
    //         "divqrsw",
    //         mesh_.facesInstance(),
    //         mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE,
    //         false
    //     )
    // ),

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
        dimensionedScalar("0", dimMoles/(dimLength*dimLength*dimTime), 0.0)
    ), 
    gsr_
    (
        IOobject
        (
            "gsr",
            runTimeSoil_.timeName(),
            meshSoil_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshSoil_,
        dimensionedScalar("0", dimTime/dimLength, 0.0)
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
        dimensionedScalar("0", dimPressure/dimPressure, 0.0)
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
        dimensionedScalar("0", dimTemperature, 0)
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
        dimensionedScalar("0", dimEnergy/(dimLength*dimLength), 0.0)
    ),  

    Sh_
    (
        IOobject
        (
            "Sh",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    Su_
    (
        IOobject
        (
            "Su",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    ),
    Sw_
    (
        IOobject
        (
            "Sw",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
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
    )

    // gv_leaf_
    // (
    //     IOobject
    //     (
    //         "gv_leaf",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
    // ),      
    // LAI_
    // (
    //     IOobject
    //     (
    //         "LAI",
    //         runTime_.timeName(),//"0",//runTime_.timeName(),
    //         mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_
    // ),
    // qlat_leaf_
    // (
    //     IOobject
    //     (
    //         "qlat_leaf",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
    // ),
    // qsen_leaf_
    // (
    //     IOobject
    //     (
    //         "qsen_leaf",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
    // ),
    // ga_
    // (
    //     IOobject
    //     (
    //         "ga",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    // ),
    // gs_
    // (
    //     IOobject
    //     (
    //         "gs",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    // ),
    // qrad_leaf_
    // (
    //     IOobject
    //     (
    //         "qrad_leaf",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
    // ),  
    // Tl_
    // (
    //     IOobject
    //     (
    //         "Tl",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(0,0,0,1,0,0,0), 0.0)
    // ),    
    // h_ch_
    // (
    //     IOobject
    //     (
    //         "h_ch",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(1,0,-3,-1,0,0,0), 0.0)
    // ),
    // h_cm_
    // (
    //     IOobject
    //     (
    //         "h_cm",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    // ),
    {
        Info << nl << "         Defined custom vegetation model: soil-vegetation foam" << endl;
        
        Info << nl << "Vegetation : [Status]      :: Vegetation properties : \n"
             << vegetationProperties_ << nl << endl;

        
        // Info << "Std. pressure Pstd = " << constant::standard::Pstd << endl;
        // Info << "Std. temperature Tstd = " << constant::standard::Tstd << endl;
        // Info << "Avogadro constant Na = " << constant::physicoChemical::NA << endl;
        // Info << "Planck constant h = " << constant::universal::h << endl;
        // Info << "Speed of light in vacuum c = " << constant::universal::c << endl;
        // Info << "Universal gas constant R = " << constant::physicoChemical::R << endl;
        // Info << "Stefan-Boltzmann constant sigma = " << constant::physicoChemical::sigma << endl;


        Info << nl << "Vegetation : [Status]      :: Fixed constants :"  << nl
                   << "                           :: " << rhow_ << nl
                   << "                           :: " << Ra_ << nl
                   << "                           :: " << Rv_ << nl
                   << "                           :: " << cpa_ << nl
                   << "                           :: " << Mco2_ << nl
                   << "                           :: " << Lv_ << nl
                   << "                           :: " << Mw_ << nl <<endl;

        // Check units
        Info << nl << "Vegetation : [Status]      :: Input constants : " << endl;
        assertUnit(c_, dimless);
        assertUnit(cao_, dimless);
        assertUnit(g_, dimAcceleration);
        assertUnit(C_, pow(dimTime,0.5)/dimLength);
        assertUnit(l_, dimLength);
        assertUnit(Hr_, dimLength);
        assertUnit(RAI_, dimless);
        assertUnit(r_, dimLength);
        assertUnit(beta_, dimTime);
        assertUnit(betaL_, dimless/(dimPressure*dimPressure));
                
        // Output vegetation properties
        writeVegetationProperties();
        

    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// solve vegetation model
void soilVegetationModel::solve(const volVectorField& U, const rhoThermo& thermo, const volScalarField& w, const volScalarField& ws, const volScalarField& pc, const volScalarField& Kl)
{
    const volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& rho = thermo.rho();
    
    // If Air time is updated:
    // Recalculate soil properties
    // Recalculate radiation field
    if (mesh_.time().value() != internalTime)
    {
        Info << "Vegetation : [Status]      :: Updating internal clock" << endl;

        Info << "Vegetation : [Soil]        :: Calculating soil properties" << endl;
        calc_soilProperties(Kl);

        Info << "Vegetation : [MWUE]        :: Calculating Marginal water use efficiency" << endl;
        calc_marginalWUE();

        Info << "Vegetation : [Radiation]   :: Updating radiation field" << endl;
        calc_radiation();

        // Initialize leaf temperature
        forAll(Tl_, cellI)
            if (LAD_[cellI] > 10*SMALL)
                Tl_[cellI] = T[cellI];

        // Update internal clock
        internalTime = mesh_.time().value();
        Info << "Vegetation : [Status]      :: Current time: " << internalTime << endl;

    }
    else
    {
        Info << "Vegetation : [Status]      :: Soil properties and radiation already calculated." << endl;
    }

    // Magnitude of velocity
    volScalarField magU("magU", mag(U));
    
    
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            // // Initial leaf temperature
            // if (i==1)
            //     Tl_[cellI] = T[cellI];

            // Calculate vapour pressure deficit [Pa/Pa]
            VPD_[cellI] = calc_VPD(Tl_[cellI], p[cellI], w[cellI]);

            // Calculate aerodynamic resistance to CO2 [mol/m2 s]
            ga_[cellI] = calc_ga(magU[cellI], rho[cellI]);
        
        }
    }

    psi_L_.value() = -1e6;

    // Correct boundary conditions
    VPD_.correctBoundaryConditions();
    ga_.correctBoundaryConditions();

    // Export
    writeVegetationProperties();
}




// void soilVegetationModel::solve(const volVectorField& U, const rhoThermo& thermo, const volScalarField& w, const volScalarField& ws)
// {
//     const volScalarField& p = thermo.p();
//     const volScalarField& T = thermo.T();
//     const volScalarField& rho = thermo.rho();
//     double pvsat, pv, gnet;
//     double Ra = 287.042; // J/kg K
//     double Rv = 461.524;

//     // Initialize radiation and leaf temperature if internalTime is changed
//     if (mesh_.time().value() != internalTime)
//     {
//         Info << "Vegetation: [Radiation] :: updating radiation field" 
//              << ", internal time = " << internalTime << endl;
//         calc_radiation();

//         // Initialize leaf temperature
//         forAll(Tl_, cellI)
//             if (LAD_[cellI] > 10*SMALL)
//                 Tl_[cellI] = T[cellI];
//     }    

//     Info << "Vegetation: [Leaf Energy Balance] :: updating leaf temperature "
//          << ", internal time = " << internalTime << endl;

//     // Magnitude of velocity
//     volScalarField magU("magU", mag(U));

//     // Bounding velocity
//     // bound(magU, UMin_); // no longer needed

//     // solve aerodynamic, stomatal resistance
//     //resistance(U,T);
//     volScalarField new_Tl("new_Tl", Tl_);

//     // 
//     dimensionedScalar mean_ws("mean_ws", fvc::domainIntegrate(ws*RAD_)/fvc::domainIntegrate(RAD_));
//     Info << "Vegetation: [soil moisture] :: mean ws: " << mean_ws << endl;
//     Info << "Vegetation: [soil moisture] :: PWP ws: " << wPWP_ << endl;



//     scalar maxError, maxRelError;
//     int i;

//     // solve leaf temperature, iteratively.
//     int maxIter = 500;
//     for (i=1; i<=maxIter; i++)
//     {
//         // Solve aerodynamic, stomatal resistance
//         //std::clock_t startTime = std::clock();
//         //resistance(magU, T, w, new_Tl);
//         //Info << "It took " << (std::clock() - startTime ) / (double) CLOCKS_PER_SEC << " second(s)." << endl;

//         forAll(LAD_, cellI)
//         {
//             if (LAD_[cellI] > 10*SMALL)
//             {
//                 // Initial leaf temperature
//                 // if (i==1)
//                 //     Tl_[cellI] = T[cellI];//*0. + 300.;//T[cellI];

//                 // vapour pressure
//                 pv = p[cellI] * w[cellI] / (Ra/Rv + w[cellI]); // vapour pressure

//                 // saturation vapour pressure
//                 pvsat = calc_pvsat(Tl_[cellI]);

//                 // aerodynamic resistance
//                 //ra_[cellI] = calc_resistance_aerodynamic(magU[cellI]);
//                 // Aerodynamic conductance (boundary-layer conductance) m/s
//                 ga_[cellI] = calc_conductance_aerodynamic(magU[cellI]);

//                 // stomatal resistance
//                 // rs_[cellI] = calc_resistance_stomatal(pv, pvsat, T[cellI], cellI);
//                 // stomatal conductance m/s
//                 gs_[cellI] = calc_conductance_stomatal(pv,pvsat,T[cellI],cellI, mean_ws);

//                 // net heat/vapor conductance m/s
//                 gnet = (gs_[cellI] * ga_[cellI]) / (gs_[cellI] + ga_[cellI]) ;

//                 // convective heat transfer coefficient
//                 //h_ch_[cellI] = (2.0*rho[cellI]*cpa_.value())/ra_[cellI];
//                 h_ch_[cellI] = 2.0 * ga_[cellI] * rho[cellI] * cpa_.value();

//                 // convective mass transfer coefficient
//                 //h_cm_[cellI] = (rho[cellI]*Ra)/(p[cellI]*Rv*(ra_[cellI]+rs_[cellI]));
//                 h_cm_[cellI] = gnet * (rho[cellI]*Ra / (p[cellI]*Rv) ); //*(ra_[cellI]+rs_[cellI]));

//                 // Calculate transpiration rate (mass flux rate)
//                 gv_leaf_[cellI] = nEvapSides_.value() * h_cm_[cellI] * (pvsat - pv);
                
//                 // Calculate latent heat flux
//                 qlat_leaf_[cellI] = lambda_.value() * gv_leaf_[cellI];

//                 // Calculate sensible heat flux
//                 qsen_leaf_[cellI] = h_ch_[cellI] * (Tl_[cellI] - T[cellI]);

//                 // Calculate new leaf temperature
//                 //new_Tl[cellI] = T[cellI] + (qrad_leaf_[cellI] - qlat_leaf_[cellI]) * (ra_[cellI]/(2.0*rho[cellI]*cpa_.value()));
//                 new_Tl[cellI] = T[cellI] + (qrad_leaf_[cellI] - qlat_leaf_[cellI]) / h_ch_[cellI];

//             }
//         }
        
//         // Check rel. L-infinity error
//         maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
//         maxRelError = maxError/gMax(mag(new_Tl.internalField()));

//         // Info
//         Info << "Vegetation: [Leaf Energy Balance] :: Iteration = " << i
//              << " max. Tl = " << gMax(new_Tl)
//              << ", error Tl = "   << maxError
//              << " LEB error = " << gMax(qrad_leaf_.internalField() - qsen_leaf_.internalField() - qlat_leaf_.internalField())
//              << endl;
        
//         // convergence check
//         if ((maxRelError < 1e-8) && (maxError < 1e-8))
//             break;
//         else
//             forAll(Tl_, cellI)
//                 Tl_[cellI] = 0.5*Tl_[cellI]+0.5*new_Tl[cellI]; // stabilized // update leaf temp.       
//     }

//     // Correct boundary conditions
//     ga_.correctBoundaryConditions();
//     gs_.correctBoundaryConditions();
//     gv_leaf_.correctBoundaryConditions();
//     qlat_leaf_.correctBoundaryConditions();
//     qsen_leaf_.correctBoundaryConditions();

//     // Iteration info
//     Info << "Vegetation: [Leaf Energy Balance] :: Final residual = " << maxError
//          << ", Final rel. residual = " << maxRelError
//          << ", No Iterations " << i 
//          << endl;
    
//     Info << "Vegetation: [Leaf Energy Balance] :: int. a*qr = " << fvc::domainIntegrate(qrad_leaf_ * LAD_).value()
//          << ", int. a*qs = " << fvc::domainIntegrate(qsen_leaf_ * LAD_).value()
//          << ", int. a*ql = " << fvc::domainIntegrate(qlat_leaf_ * LAD_).value() 
//          << endl;

// }

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> soilVegetationModel::Sh()
{
    //Sh_ = LAD_ * qsen_leaf_;
    //Sh_.correctBoundaryConditions();
    return Sh_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> soilVegetationModel::Su(volScalarField& rho, volVectorField& U)
{
    // Calculate
    Su_ = -Cf_*rho*mag(U)*U;
    // Correct boundary conditions
    Su_.correctBoundaryConditions();

    return fvm::SuSp(-Cf_*rho*mag(U), U);
}

// return humidity source term (moisture content)
tmp<volScalarField> soilVegetationModel::Sw()
{
    //Sw_ = LAD_ * gv_leaf_;
    //Sw_.correctBoundaryConditions();
    return Sw_;
}

// return soil moisture source term
tmp<volScalarField> soilVegetationModel::Sws(volScalarField& Kl, volScalarField& Cl, volScalarField& pc)
{
    // Hydraulic liquid diffusivity, m2/s
    //volScalarField Dl("Dl", Kl/Cl);

    // Net transpiration rate kg/s
    //dimensionedScalar mtrans("mtrans", fvc::domainIntegrate(LAD_ * gv_leaf_));
 
    // source term for water uptake due to roots kg/(m3s)
    //Sws_ = - (mtrans * RAD_ * Dl) / fvc::domainIntegrate(RAD_ * Dl);
    //Spc_ = - (0.0 * mtrans * RAD_ * Dl) / fvc::domainIntegrate(RAD_ * Dl);
    
    //Info << "Vegetation: [soil root uptake] :: Sws_ max = " << gMax(Sws_) 
    //     << ", min: " << gMin(Sws_) << endl;

    return Sws_;
}


// -----------------------------------------------------------------------------

bool soilVegetationModel::read()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam



// /*
// double soilVegetationModel::calc_resistance_aerodynamic(const double& magU)
// {
//     return (C_.value() * pow(l_.value() / magU, 0.5));
// }
// */

// double soilVegetationModel::calc_rhosat(const double& T)
// {
//     return calc_pvsat(T)/(461.5*T);
// }

// // void soilVegetationModel::calc_radiation()
// // {
// //     label timestepsInADay_ = divqrsw.size();

// //     Time& time = const_cast<Time&>(mesh_.time());
// //     Info << "Vegetation: [Radiation] :: time.value(): " << time.value();

// //     label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
// //     Info << ", 1 timestep: " << timestep;
// //     timestep = timestep % timestepsInADay_;
// //     Info << ", 2 timestep: " << timestep << endl;

// //     scalarList divqrswi =  divqrsw[timestep];

// //     // radiation density inside vegetation
// //     forAll(Cf_, cellI)
// //         if (Cf_[cellI] > 10*SMALL)
// //             qrad_leaf_[cellI] = (- divqrswi[cellI] )/ LAD_[cellI];

// //     qrad_leaf_.correctBoundaryConditions();

// // }





// double soilVegetationModel::calc_conductance_stomatal(const double& pv, const double& pvsat, const double& T, const int& cellI, const dimensionedScalar& mean_ws)
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

//     if (isDayTime_)
//     {
//         // constant Stomatal resistance
//         // case STOMATAL_RESISTANCE_CONSTANT: return rsMin_.value();
//         double rs = rsMin_.value()*(0.971 + pow(wPWP_.value()/mean_ws.value(), 2));
//         //return 1.0/rsMin_.value();    
//         return 1.0/rs;    
//     }
//     else
//     {
//         return SMALL;
//     }
    
//     // Bruse and Fleer (1998)
//     // return rsMin_.value()* // (R_sw_max / (0.03*R_sw_max + R_sw) + f_grow + pow(eta_wilt/eta, 2));
    
// }



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
