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

Class
    Foam::soilVegetationModel

Description
    soil Vegetation model

Author
    Lento Manickathan mlento@ethz.ch

SourceFiles
    vegetationModel.C

\*---------------------------------------------------------------------------*/

#ifndef soilVegetationModel_H
#define soilVegetationModel_H

#include "fvc.H"
#include "rhoThermo.H"
#include "fvm.H"
#include "IOdictionary.H"
#include "bound.H"
#include "scalarListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"
#include <ctime>
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class vegetationModel Declaration
\*---------------------------------------------------------------------------*/

class soilVegetationModel
:
    public IOdictionary
{

private:

    // Private data
    
        // disallow copy construct
        soilVegetationModel(const soilVegetationModel&);

        // disallow default bitwise assignment
        void operator=(const soilVegetationModel&);

protected:

    // Protected data


        // Model property
    
            dictionary vegetationProperties_;
    
    
            //- Reference to air time class
            const Time& runTime_;

            //- Reference to soil time class       
            const Time& runTimeSoil_;

            //- Reference to air mesh
            const fvMesh& mesh_;

            //- Reference to soil mesh
            const fvMesh& meshSoil_;

            //- Internal clock
            scalar internalTime;

            //- Minimum velocity threshold
            const dimensionedScalar UMin_;

            //- Min. threshold for LAD, RAD
            const scalar minThreshold;
        
            //bool isDayTime_;            // check if day time

        // Model coefficients
    
            // UNITS = [kg m s K mol A cd]

            //- Density of water [kg/m3]
            const dimensionedScalar rhow_;

            //- Standard pressure [Pa, kg/m s2]
            const dimensionedScalar Pstd_;

            //- Gas constant of dry air [J/kg K]
            const dimensionedScalar Ra_;

            //- Gas constant of water vapour [J/kg K]
            const dimensionedScalar Rv_;

            //- Specific heat capacity of air at const. pressure [J/kg K]
            const dimensionedScalar cpa_;

            //- Molar mass of water [kg/mol]
            const dimensionedScalar Mw_;

            //- Molar mass of CO2 [kg/mol]
            const dimensionedScalar Mco2_;

            //- Latent heat of vaporization of water
            const dimensionedScalar Lv_;

            //- Ambient CO2 Concentration [mol/mol]
            const dimensionedScalar c_;

            //- Ambient O2 Concentration [mol/mol]
            const dimensionedScalar cao_;

            //- Gravitational acceleration [m/s2]
            const dimensionedVector g_;

            //- Leaf boundary layer conductance proportionality factor
            const dimensionedScalar C_;

            //- Characteristic size of plant leaf [m]
            const dimensionedScalar l_;

            //- Root depth [m]
            const dimensionedScalar Hr_;

            //- Root area index [m2/m2]
            const dimensionedScalar RAI_;

            //- Root radius [m]
            const dimensionedScalar r_;

            //- Soil conductance [s]
            const dimensionedScalar beta_;

            //- Initial 24h average leaf water potential [Pa]
            const dimensionedScalar psi_L24_0_;

            //- Timesteps in a day
            const label timestepsInADay_;

            //- Maximum marginal water use efficiency (mol/mol)
            const dimensionedScalar lambda_max_;

            //- CO2 at maximum. marginal WUE (mol/mol)
            const dimensionedScalar c_star_;

            //- Marginal water use efficiency Coefficient
            const dimensionedScalar betaL_;

            //- Leaf water potential at maximum WUE
            const dimensionedScalar psi_Lmax_;

            //- Leaf water potential [Pa]
            dimensionedScalar psi_L_;

            //- Root water potential [Pa]
            dimensionedScalar psi_R_;

            //- 24h average Leaf water potential [Pa]
            scalarList psi_L24_;

            //- Marginal water use efficiency
            dimensionedScalar lambda_;

            //- Varying vegetation properties
            IOdictionary varyingVegetationProperties_;

            //dimensionedScalar a1_;          // stomatal resistance specific model parameters
            //dimensionedScalar a2_;          // stomatal resistance specific model parameters
            //dimensionedScalar a3_;          // stomatal resistance specific model parameters
            
            // //- Ambient O2 Concentration [mol/mol]
            // const dimensionedScalar cao_("cao", dimensionSet(0,0,0,0,0,0,0), 210); // 210 mmol/mol - 0.21 mol/mol
            //     dimensionedScalar("0", dimensionSet(1,-2,-1,0,0,0,0), 0.0)

            
            // //- Ambient CO2 Concentration [mumol/mol]
            // dimensionedScalar ca_; // 380 mumol/mol - 0.00038 mol/mol

            // //- Specific heat capacity of air at constant pressure [J/kg K]
            // dimensionedScalar cpa_;

            // //- Density of water // [kg/m3]
            // dimensionedScalar rhow_;

            // //- Gas constant air =
            // //dimensionedScalar Ra_("Ra", [0, 0, 0 ,0 0 0 0 ], 212.1);
            
            
            // //dimensionedScalar D0_;          // vapour pressure corresponding to minimal resistance
            
            // dimensionedScalar nEvapSides_;  // number of sides, leaf evaporates from
            
            // //dimensionedScalar H_;           // Height of vegetation
            
            // //dimensionedScalar kc_;          // extinction coefficient of radiation
            
            
            // //dimensionedScalar Rg0_;         // global radiation at top of canopy
            // //dimensionedScalar Rl0_;         // long-wave radiation at top of canopy
            
            // //dimensionedScalar rhoa_;        // denisty of air
            // dimensionedScalar rsMin_;       // minimum stomatal resistance
            
            //dimensionedScalar TlMin_;       // lower limit of T

            //- lower limit of U, for stability
            

            // dimensionedScalar lambda_;      // latent heat of vaporization
            // //dimensionedScalar wPWP_;      // permanent wilting point moisture content
            
            // dimensionedScalar mtrans_;      // net transpirated mass of water vapour

            const scalarListIOList divqrsw;
    

        // Fields

            //- Vegetation net drag = cd*LAD [m2/m3]
            volScalarField Cf_;                 

            //- Leaf area density [m2/m3]
            volScalarField LAD_;        

            //- Aerodynamic conductance to CO2 [mol/m2 s]
            volScalarField ga_;
    
            //- Effective soil-root conductance to water [s/m]
            volScalarField gsr_;

            //- Leaf to air vapor pressure deficit [Pa/Pa]
            volScalarField VPD_;             

            //- Root area density [m2/m3]
            volScalarField RAD_;                

            //- Leaf temperature [K]
            volScalarField Tl_;                

            //- Short + long wave radiation absorbed by leaves (W/m2)
            volScalarField qrad_leaf_;         
            
            //- Energy source term [kg/m s3]
            volScalarField Sh_;        

            //- Momentum source term [kg/m2 s2]
            volVectorField Su_;

            //- Absolute mositure source term [kg/m3 s]
            volScalarField Sw_;                 
            
             //- Soil moisture source term [kg/m3 s]
            volScalarField Sws_;

            // //volScalarField E_;                // transpiration rate
            // volScalarField gv_leaf_;            // transpiration rate
            // //volScalarField pv_;               // water vapor pressure
            // //volScalarField pvsat_;            // saturated water vapor pressure
            
            

            // volScalarField LAI_;                // leaf area index
            // //volScalarField wsat_;             // saturated specific humidity
            // //volScalarField Qlat_;             // latent heat flux
            //volScalarField qlat_leaf_;          // latent heat flux
            // //volScalarField Qsen_;             // sensible heat flux
            // volScalarField qsen_leaf_;          // sensible heat flux
            
            // volScalarField gs_;                 // stomatal resistance
            // //volScalarField rhosat_;           // saturated density field
            // //volVectorField Rg_;               // global radiation
            // //volScalarField Rn_;               // global radiation density in volume
            
            // volScalarField Tl_;                 // temperature of leaf

            // volScalarField h_ch_;               // convective heat transfer coefficient
            // volScalarField h_cm_;               // convective mass transfer coefficient


    // Protected Member Functions

        //- Assert the units are as required
        template<class tmpClass>
        void assertUnit(const tmpClass& sourceVar, const dimensionSet& targetUnit);
        //virtual void assertUnit(const dimensionedScalar& sourceVar, const dimensionSet& targetUnit);
        
        //- Write time-varying vegetation properties
        virtual void writeVegetationProperties();

        // //- Return saturated density of water vapour
        // virtual double calc_rhosat(const double& T);

        //- Return saturated water vapour pressure [Pa]
        virtual double calc_pvsat(const double& T);

        //- Return vapour pressure [Pa]
        virtual double calc_pv(const double& p, const double& w);

        //- Return vapour pressure deficit [Pa/Pa]
        virtual double calc_VPD(const double& T, const double& p, const double& w);

        //- Return stomatal condutance to CO2 (mol/m2 s)
        virtual double calc_ga(const double& magU, const double& rho);

        //- Calculate soil properties
        virtual void calc_soilProperties(const volScalarField& Kl);


        //- Calculate marginal water use efficiency
        virtual void calc_marginalWUE();

        //- calculate radiation
        virtual void calc_radiation();
        
        /*
        //- calculate aerodynamic, stomatal resistances
        //virtual void resistance(const volScalarField& magU, const volScalarField& T, const volScalarField& w, const volScalarField& Tl);

        // solve aerodynamic resistance
        //virtual double calc_resistance_aerodynamic(const double& magU);
        virtual double calc_conductance_aerodynamic(const double& magU);

        // solve stomtal resistance
        //virtual double calc_resistance_stomatal(const double& pv, const double& pvsat, const double& T, const int& cellI);
        //virtual double calc_conductance_stomatal(const double& pv, const double& pvsat, const double& T, const int& cellI);
        virtual double calc_conductance_stomatal(const double& pv, const double& pvsat, const double& T, const int& cellI, const dimensionedScalar& mean_ws);
        */

public:
    
    //-Runtime type information
    TypeName("soilVegetationModel");

    // Constructors
    
        //- Construct from components
        soilVegetationModel
        (
            const volVectorField& U,    //- Fluid velocity [m/s]
            const rhoThermo& thermo,    //- Fluid thermodynamic properties
            const volScalarField& w,    //- Fluid absolute humidity [kg/kg]
			const volScalarField& Ts,   //- Soil temperature [K]
            const volScalarField& ws,   //- Soil moisture [kg/m3]
            const volScalarField& pc    //- Soil capillary pressure [Pa, kg/ms2]
        );

    //- Destructor
    virtual ~soilVegetationModel()
    {}

    // Member Functions
        
        //- Return vegetation properties dictionary
        dictionary vegetationProperties() const
        {
            return vegetationProperties_;
        }
    
        // //- Solve leaf energy balance
        //virtual void solve(volVectorField&U, volScalarField& T, volScalarField& w);//, volScalarField& Tl_);
        //virtual void solve(const volVectorField&U, const rhoThermo& thermo, const volScalarField& w);//, volScalarField& Tl_);
        virtual void solve(const volVectorField&U, const rhoThermo& thermo, const volScalarField& w, const volScalarField& ws, const volScalarField& pc, const volScalarField& Kl);

        // Source terms

        //- Return vegetation energy source
        virtual tmp<volScalarField> Sh();

        //- Return vegetation specific humidity source
        virtual tmp<volScalarField> Sw();

        //- Return vegetation momentum source
        virtual tmp<fvVectorMatrix> Su(volScalarField& rho, volVectorField& U);

        //- Return soil moisture source term
        virtual tmp<volScalarField> Sws(volScalarField& Kl, volScalarField& Cl, volScalarField& pc); 

        //- Read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif


// ************************************************************************* //
