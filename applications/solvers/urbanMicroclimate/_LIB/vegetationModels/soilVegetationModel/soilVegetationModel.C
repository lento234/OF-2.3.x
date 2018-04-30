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
#include "leaffluxes.H"
#include "moisture.H"
#include "plantfluxes.H"
#include "radiation.H"

using namespace Foam::constant;
//using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(soilVegetationModel, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Assert the units are as required
template<class tmpClass>
void soilVegetationModel::assertDimensions(const tmpClass& sourceVar, const dimensionSet& targetUnit)
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

//- Write time-varying vegetation properties
void soilVegetationModel::writeVegetationProperties()
{
    // Update dictionary
    varyingVegetationProperties_.set("psi_L", psi_L_);
    varyingVegetationProperties_.set("psi_R", psi_R_);
    varyingVegetationProperties_.set("lambda_", lambda_);

    //Info
    Info << "Vegetation : [Status]      :: Writing varying vegetation properties" << endl;

    // Export to file
    varyingVegetationProperties_.regIOobject::write();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

soilVegetationModel::soilVegetationModel
(
    const volVectorField& U,    //- Fluid velocity [m/s]
    const rhoThermo& thermo,    //- Fluid thermodynamic properties
    const volScalarField& w,    //- Fluid absolute humidity [kg/kg]
    const volScalarField& c,    //- Fluid CO2 concentration [mol/mol]
    //const volScalarField& Ts,   //- Soil temperature [K]
    //const volScalarField& ws,   //- Soil moisture [kg/m3]
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
    runTimeSoil_(pc.time()),
    mesh_(U.mesh()),
    meshSoil_(pc.mesh()),
    UMin_("UMin", dimVelocity, SMALL),
    minThreshold(10*SMALL),

    rhow_("rhow", dimDensity, 1000.0),
    Pstd_("Pstd", dimPressure, 101300),
    Ra_("Ra", dimGasConstant, 287.042),
    Rv_("Rv", dimGasConstant, 461.524),
    cpa_("cpa", dimSpecificHeatCapacity, 1003.5),
    Mw_("Mw", dimMass/dimMoles, 0.01802),
    Mco2_("Mco2", dimMass/dimMoles, 0.04401),
    Mair_("Mair", dimMass/dimMoles, 0.02897),
    Lv_("Lv", dimEnergy/dimMass, 2.5e6),
    a_(1.6),
    E_PAR_("E_PAR", dimEnergy/dimMoles, 2.24e5),


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

    ca_
    (
        vegetationProperties_.subDict("environmentalCoeffs").lookup("ca")
    ),
    cao_
    (
        vegetationProperties_.subDict("environmentalCoeffs").lookup("cao")
    ),
    g_
    (
        vegetationProperties_.subDict("environmentalCoeffs").lookup("g")
    ),
    C_
    (
        vegetationProperties_.subDict("leafCoeffs").lookup("C")
    ),
    l_
    (
        vegetationProperties_.subDict("leafCoeffs").lookup("l")
    ),
    Hr_
    (
        vegetationProperties_.subDict("rootCoeffs").lookup("Hr")
    ),
    RAI_
    (
        vegetationProperties_.subDict("rootCoeffs").lookup("RAI")
    ),
    r_
    (
        vegetationProperties_.subDict("rootCoeffs").lookup("r")
    ),    
    beta_
    (
        vegetationProperties_.subDict("rootCoeffs").lookup("beta")
    ),
    psi_L24_0_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("psi_L24_0")
    ),
    timestepsInADay_
    (
        vegetationProperties_.lookupOrDefault("timestepsInADay", 24)
    ),
    lambda_max_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("lambda_max")
    ),    
    ca_star_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("ca_star")
    ),    
    betaL_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("betaL")
    ),  
    psi_Lmax_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("psi_Lmax")
    ),
    gx_max_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("gx_max")
    ),
    d_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("d")
    ),
    cx_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookupOrDefault("cx", 2.0)
    ),
    Ax_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("Ax")
    ),
    s_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookupOrDefault("s", 0.7)
    ),
    gsn_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("gsn")
    ),
    Vcmax25_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("Vcmax25")
    ),
    Kc25_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("Kc25")
    ),
    Ko25_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("Ko25")
    ),
    gammac_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookupOrDefault("gammac", 0.074)
    ),
    gammao_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookupOrDefault("gammao", 0.018)
    ),
    kco_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("kco")
    ),
    ko_
    (
        vegetationProperties_.subDict("stomataCoeffs").lookup("ko")
    ),        
    rPAR_
    (
        vegetationProperties_.subDict("radCoeffs").lookupOrDefault("rPAR", 0.5)
    ),
    gammaPAR_
    (
        vegetationProperties_.subDict("radCoeffs").lookupOrDefault("gammaPAR", 0.015)
    ),    
    nEvapSides_
    (
        vegetationProperties_.subDict("leafCoeffs").lookupOrDefault("nEvapSides", 1.0)
    ),


    gabs_("g", mag(g_)),
    alpha_("alpha", pow(Hr_/RAI_, 0.5)/pow(2.0*r_, 0.5)),


    psi_L_("psi_L", dimPressure, 0.0),
    psi_R_("psi_R", dimPressure, 0.0),
    psi_L24_(timestepsInADay_,psi_L24_0_.value()),
    lambda_("lambda", dimMoles/dimMoles, 0.0),
    internalTime(-1.0),


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
    qrsw_leaf_
    (
        IOobject
        (
            "qrsw_leaf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimEnergy/(dimLength*dimLength), 0.0)
    ),  
    qrlw_leaf_
    (
        IOobject
        (
            "qrlw_leaf",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimEnergy/(dimLength*dimLength), 0.0)
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
        dimensionedScalar("0", dimMoles/(dimLength*dimLength*dimTime), 0.0)
    ),  
    gc_eff_
    (
        IOobject
        (
            "gc_eff",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimMoles/(dimLength*dimLength*dimTime), 0.0)
    ),  
    gv_eff_
    (
        IOobject
        (
            "gv_eff",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimMoles/(dimLength*dimLength*dimTime), 0.0)
    ),      
    An_
    (
        IOobject
        (
            "An",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimMoles/(dimLength*dimLength*dimTime), 0.0)
    ), 
    ci_
    (
        IOobject
        (
            "ci",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimMoles/dimMoles, 0.0)
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
        dimensionedScalar("0", dimMass/(dimLength*dimLength*dimTime), 0.0)
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
        dimensionedScalar("0", dimPower/(dimLength*dimLength), 0.0)
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
        dimensionedScalar("0", dimPower/(dimLength*dimLength), 0.0)
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
        dimensionedScalar("0", dimPower/(dimLength*dimLength*dimLength), 0.0)
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
    Sc_
    (
        IOobject
        (
            "Sc",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", (dimMass/(dimLength*dimLength*dimLength*dimTime))*dimMoles/dimMoles, 0.0)
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
        assertDimensions(ca_, dimless);
        assertDimensions(cao_, dimless);
        assertDimensions(g_, dimAcceleration);
        assertDimensions(C_, pow(dimTime,0.5)/dimLength);
        assertDimensions(l_, dimLength);
        assertDimensions(Hr_, dimLength);
        assertDimensions(RAI_, dimless);
        assertDimensions(r_, dimLength);
        assertDimensions(beta_, dimTime);
        assertDimensions(betaL_, dimless/(dimPressure*dimPressure));
        assertDimensions(psi_Lmax_, dimPressure);
        assertDimensions(gx_max_, dimless/dimTime);
        assertDimensions(d_, dimPressure);
        assertDimensions(Ax_, dimLength*dimLength);
        assertDimensions(gsn_, dimMoles/(dimLength*dimLength*dimTime));        
        assertDimensions(Vcmax25_, dimMoles/(dimLength*dimLength*dimTime));        
        assertDimensions(Kc25_, dimMoles/dimMoles);
        assertDimensions(Ko25_, dimMoles/dimMoles);
        assertDimensions(kco_, dimless/dimTime);
        assertDimensions(ko_, dimless/dimTime);

        // read relaxation factor for Tl - aytac
        //dictionary relaxationDict = mesh_.solutionDict().subDict("relaxationFactors");
        //scalar Tl_relax = relaxationDict.lookupOrDefault<scalar>("Tl", 0.5);

        

        // Output vegetation properties
        writeVegetationProperties();

    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// solve leaf energy balance 
void soilVegetationModel::solve(const volVectorField& U, const rhoThermo& thermo, const volScalarField& w, const volScalarField& c)
{
    const volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& rho = thermo.rho();
    
    // If Air time is updated: Recalculate soil properties + Recalculate radiation field
    if (mesh_.time().value() != internalTime)
    {
        Info << "Vegetation : [Status]      :: Updating internal clock" << endl;

        Info << "Vegetation : [MWUE]        :: Calculating Marginal water use efficiency" << endl;
        calc_marginalWUE();

        Info << "Vegetation : [Radiation]   :: Updating radiation field" << endl;
        calc_radiation();

        // Initialize leaf temperature
        forAll(Tl_, cellI)
            if (LAD_[cellI] > minThreshold)
                Tl_[cellI] = T[cellI];

        // Update internal clock
        internalTime = mesh_.time().value();
        Info << "Vegetation : [Status]      :: Current time: " << internalTime << endl;

    }
    else
    {
        Info << "Vegetation : [Status]      :: MWUE and Radiation already calculated." << endl;
    }

    
    // New guess of leaf temperature
    volScalarField new_Tl("new_Tl", Tl_);
    
    scalar maxError, maxRelError;
    label i;

    //- Calculate aerodynamic resistance (mol / m2 s)
    calc_ga(U, rho);

    scalar CHTC, CMTC;
    // solve leaf temperature, iteratively.
    label maxIter = 500;
    for (i=1; i<=maxIter; i++)
    {

        // Loop through all the leaf cells
        forAll(LAD_, cellI)
        {
            // Only leaf cells
            if (LAD_[cellI] > minThreshold)
            {
                // Calculate vapour pressure deficit [Pa/Pa]
                VPD_[cellI] = calc_VPD(Tl_[cellI], p[cellI], w[cellI]);

                // Solve assimilation (An, gs, gc_eff, gv_ef, ci)
                solve_assimilation(Tl_[cellI], VPD_[cellI], c[cellI], cellI);

                // convective heat transfer coefficient (W / m2 K)
                CHTC = 2.0 * cpa_.value() * Mco2_.value() * ga_[cellI]; 

                // convective mass transfer coefficient (kg / m2 s)
                CMTC = Mw_.value() * gv_eff_[cellI];

                // Calculate transpiration rate (mass flux rate)
                gv_leaf_[cellI] = nEvapSides_ * CMTC * VPD_[cellI];
                
                // Calculate latent heat flux (W / m2)
                qlat_leaf_[cellI] = Lv_.value() * gv_leaf_[cellI];

                // Calculate sensible heat flux (W / m2)
                qsen_leaf_[cellI] = CHTC * (Tl_[cellI] - T[cellI]);

                // Calculate new leaf temperature
                new_Tl[cellI] = T[cellI] + (qrad_leaf_[cellI] - qlat_leaf_[cellI]) / CHTC;
            }
        }
        
        // Check rel. L-infinity error
        maxError = gMax(mag(new_Tl.internalField()-Tl_.internalField()));
        maxRelError = maxError/gMax(mag(new_Tl.internalField()));

        // Info
        // Info << "Vegetation: [Leaf Energy Balance] :: Iteration = " << i
        //      << " max. Tl = " << gMax(new_Tl)
        //      << ", error Tl = "   << maxError
        //      //<< " LEB error = " << gMax(qrad_leaf_.internalField() - qsen_leaf_.internalField() - qlat_leaf_.internalField())
        //      << endl;
        
        // convergence check
        if ((maxRelError < 1e-12) && (maxError < 1e-12))
            break;
        else
            forAll(Tl_, cellI)
                Tl_[cellI] = 0.5*Tl_[cellI]+0.5*new_Tl[cellI]; // stabilized // update leaf temp.       
    }

    // Correct boundary conditions
    VPD_.correctBoundaryConditions();
    gs_.correctBoundaryConditions();
    gc_eff_.correctBoundaryConditions();
    gv_eff_.correctBoundaryConditions();
    ci_.correctBoundaryConditions();
    An_.correctBoundaryConditions();
    gv_leaf_.correctBoundaryConditions();
    qlat_leaf_.correctBoundaryConditions();
    qsen_leaf_.correctBoundaryConditions();

    // Iteration info
    Info << "Vegetation : [LEB]         :: Final residual = " << maxError
         << ", Final rel. residual = " << maxRelError
         << ", No Iterations " << i 
         << endl;
    
    Info << "Vegetation : [LEB]         :: int. a*qr = " << fvc::domainIntegrate(qrad_leaf_ * LAD_).value()
         << ", int. a*qs = " << fvc::domainIntegrate(qsen_leaf_ * LAD_).value()
         << ", int. a*ql = " << fvc::domainIntegrate(qlat_leaf_ * LAD_).value() 
         << endl;


    // Export time-dependent vegetation properties
    writeVegetationProperties();
}

// -----------------------------------------------------------------------------

// return energy source term
tmp<volScalarField> soilVegetationModel::Sh()
{
    Sh_ = LAD_ * qsen_leaf_;
    
    // Correct boundary condition
    Sh_.correctBoundaryConditions();

    return Sh_;
}

// solve & return momentum source term (explicit)
tmp<fvVectorMatrix> soilVegetationModel::Su(volScalarField& rho, volVectorField& U)
{
    // Calculate
    Su_ = -Cf_*rho*mag(U)*U;
    
    // Correct boundary condition
    Su_.correctBoundaryConditions();

    return fvm::SuSp(-Cf_*rho*mag(U), U);
}

// return humidity source term (moisture content)
tmp<volScalarField> soilVegetationModel::Sw()
{

    Sw_ = LAD_ * gv_leaf_;

    //Correct boundary condition
    Sw_.correctBoundaryConditions();
    return Sw_;
}

// return co2 source term  (mol/mol * kg/m3 s)
tmp<volScalarField> soilVegetationModel::Sc()
{
    Sc_ = (LAD_ * An_ * Mair_);

    //Correct boundary condition
    Sc_.correctBoundaryConditions();
    return Sc_;
}

// return soil moisture source term (WIP)
tmp<volScalarField> soilVegetationModel::Sws(volScalarField& Kl, volScalarField& pc)
{
    // Hydraulic liquid diffusivity, m2/s
    //volScalarField Dl("Dl", Kl/Cl);

    /*
    // Solve Soil-Plant-Atmosphere Continuum
    solve_SPAC(const volScalarField& pc, const volScalarField& Kl)

    
    // source term for water uptake due to roots kg/(m3s)
    Sws_ = - gsr_ * (psi_S - psi_R_) * RAD_;
    
    Sws_.correctBoundaryConditions();
    */

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