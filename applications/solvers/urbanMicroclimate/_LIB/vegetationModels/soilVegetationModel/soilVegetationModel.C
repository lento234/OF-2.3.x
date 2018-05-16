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
#include "soilroot.H"

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
    // else //debug
    // {
    //     Info << sourceVar << endl;
    // }
}

//- Write time-varying vegetation properties
void soilVegetationModel::writeVegetationProperties()
{
    
    //Info
    Info << "Vegetation : [Status]      :: Writing varying vegetation properties" << endl;

    // Update dictionary

    // Scalar properties
    varyingVegetationProperties_.set("psi_L", psi_L_); // leaf water potential
    varyingVegetationProperties_.set("psi_R", psi_R_); // root water potential
    varyingVegetationProperties_.set("psi_L24av", dimensionedScalar("psi_L24av", dimPressure, average(psi_L24_)) ); // water potential 24 hr
    varyingVegetationProperties_.set("lambda", lambda_); // WUE
    varyingVegetationProperties_.set("Qp", Qp_); // PAR flux density
    varyingVegetationProperties_.set("E", dimensionedScalar("E", E_)); // Net vapour flux
    varyingVegetationProperties_.set("gx", gx_); // xylem conductance

    // Min, Maximum leaf temperature
    scalar Tlmax = 0.0;
    scalar Tlmin = 1000.0;
    forAll(LAD_,cellI)
    {
        if (LAD_[cellI] > minThreshold)
        {
            Tlmax = max(Tlmax, Tl_[cellI]);
            Tlmin = min(Tlmin, Tl_[cellI]);
        }
    }
    List<scalar> Tlmax_(Pstream::nProcs());
    List<scalar> Tlmin_(Pstream::nProcs());
    Tlmax_[Pstream::myProcNo()] = Tlmax;
    Tlmin_[Pstream::myProcNo()] = Tlmin;
    Pstream::gatherList(Tlmax_);
    Pstream::scatterList(Tlmax_);
    Pstream::gatherList(Tlmin_);
    Pstream::scatterList(Tlmin_);
    Tlmax = gMax(Tlmax_);
    Tlmin = gMin(Tlmin_);

    varyingVegetationProperties_.set("Tlmax", dimensionedScalar("Tlmax", dimTemperature, Tlmax)); // maxmimum leaf temperature
    varyingVegetationProperties_.set("Tlmin", dimensionedScalar("Tlmin", dimTemperature, Tlmin)); // minimum leaf temperature

    // Energy balance
    dimensionedScalar Qrad("Qrad", fvc::domainIntegrate(qrad_leaf_ * LAD_));
    dimensionedScalar Qlat("Qlat", fvc::domainIntegrate(qlat_leaf_ * LAD_));
    dimensionedScalar Qsen("Qsen", fvc::domainIntegrate(qsen_leaf_ * LAD_));
    
    varyingVegetationProperties_.set("Qrad", Qrad);
    varyingVegetationProperties_.set("Qlat", Qlat);
    varyingVegetationProperties_.set("Qsen", Qsen);

    // Assimilation Total CO2 extraction
    dimensionedScalar An("An", fvc::domainIntegrate(An_ * LAD_));
    varyingVegetationProperties_.set("An", An);

    
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
        vegetationProperties_.subDict("radCoeffs").lookupOrDefault("rPAR", 0.8)
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


    psi_L_("psi_L", dimPressure, psi_L24_0_.value()),
    psi_R_("psi_R", dimPressure, 0.0),
    psi_L24_(timestepsInADay_, psi_L24_0_.value()),
    lambda_("lambda", dimMoles/dimMoles, 0.0),
    Qp_("Qp", dimMoles/(dimLength*dimLength*dimTime), 0.0),    
    E_("E", dimMass/dimTime, 0.0),   
    gx_("gx", dimTime/dimLength, 0.0),   
    H_("H", dimLength, 0.0),
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
    gv_root_
    (
        IOobject
        (
            "gv_root",
            runTimeSoil_.timeName(),
            meshSoil_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshSoil_,
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
        
        Info << nl << nl << "Vegetation : [Status]      :: Vegetation properties : \n"
             << vegetationProperties_ << nl << endl;

        // Info << "Std. pressure Pstd = " << constant::standard::Pstd << endl;
        // Info << "Std. temperature Tstd = " << constant::standard::Tstd << endl;
        // Info << "Avogadro constant Na = " << constant::physicoChemical::NA << endl;
        // Info << "Planck constant h = " << constant::universal::h << endl;
        // Info << "Speed of light in vacuum c = " << constant::universal::c << endl;
        // Info << "Universal gas constant R = " << constant::physicoChemical::R << endl;
        // Info << "Stefan-Boltzmann constant sigma = " << constant::physicoChemical::sigma << endl;

        // Check units
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

        
        // Determine tree height
        vector nz = -(g_/gabs_).value();
        forAll(LAD_,cellI)
            if (LAD_[cellI] > minThreshold)
                H_.value() = max(mesh_.C()[cellI] & nz, H_.value());
        
        List<scalar> Hlist_(Pstream::nProcs());
        Hlist_[Pstream::myProcNo()] = H_.value();
        Pstream::gatherList(Hlist_);
        Pstream::scatterList(Hlist_);
        H_.value() = gMax(Hlist_);

        Info << "Tree height H = " << H_ << endl;        

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
    label maxIter = 500; //#TODO# : should in solution dict
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
        
        // convergence check
        if ((maxRelError < 1e-12) && (maxError < 1e-12)) //#TODO# : should in solution dict
            break;
        else
            forAll(Tl_, cellI)
                Tl_[cellI] = 0.5*Tl_[cellI]+0.5*new_Tl[cellI];// //#TODO# : should in solution dict, stabilized // update leaf temp.       
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

    //- Store net transpiration rate (kg/s) 
    E_ = fvc::domainIntegrate(LAD_ * gv_leaf_);

    // Iteration info
    Info << "Vegetation : [LEB]         :: Final residual = " << maxError
         << ", Final rel. residual = " << maxRelError
         << ", No Iterations " << i 
         << endl;
    
    Info << "Vegetation : [LEB]         :: int. a*qr = " << fvc::domainIntegrate(qrad_leaf_ * LAD_).value()
         << ", int. a*qs = " << fvc::domainIntegrate(qsen_leaf_ * LAD_).value()
         << ", int. a*ql = " << fvc::domainIntegrate(qlat_leaf_ * LAD_).value() 
         << endl;

    Info << "Vegetation : [LEB]         :: Net transpiration E = " << E_.value() << endl;
         
         

    // Export time-dependent vegetation properties
    // writeVegetationProperties();
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
    uniformDimensionedVectorField g = db().lookupObject<uniformDimensionedVectorField>("g");
    //- Calculate soil water potential
    volScalarField psi_S("psi_S", pc - ( (rhow_ * g) & meshSoil_.C() ) );
    psi_S.correctBoundaryConditions();

    // Solve Soil-Plant-Atmosphere Continuum, psi_L, psi_R, gsr
    solve_SPAC(Kl, psi_S);

    // source term for water uptake due to roots kg/(m3s)
    Sws_ = RAD_ * gv_root_;
    
    Sws_.correctBoundaryConditions();

    return Sws_;
}


// -----------------------------------------------------------------------------

bool soilVegetationModel::read()
{
    return true;
}


void soilVegetationModel::write()
{
    writeVegetationProperties();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam
