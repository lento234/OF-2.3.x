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

    // Model coefficients
    
            //- Minimum velocity threshold
            const dimensionedScalar UMin_;

            //- Min. threshold for LAD, RAD
            const scalar minThreshold;

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

            //- Molar mass of dry air [kg/mol]
            const dimensionedScalar Mair_;

            //- Latent heat of vaporization of water
            const dimensionedScalar Lv_;

            //- Relative diffusion of H2O to CO2
            const scalar a_;

            //- Energy density inside PAR (J/mol)
            const dimensionedScalar E_PAR_;


    // Model constant properties

            //- Sink term of short-wave radiation (W/m3)
            const scalarListIOList divqrsw;


            //- Ambient CO2 Concentration [mol/mol]
            const dimensionedScalar ca_;

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
            const dimensionedScalar ca_star_;

            //- Marginal water use efficiency Coefficient
            const dimensionedScalar betaL_;

            //- Leaf water potential at maximum WUE
            const dimensionedScalar psi_Lmax_;

            //- Maximum xylem condutance (1/s)
            const dimensionedScalar gx_max_;

            //- Xylem condutance constant [Pa]
            const dimensionedScalar d_;

            //- Xylem condutance constant 
            const scalar cx_;

            //- Xylem area [m2]
            const dimensionedScalar Ax_;

            //- Stomatal condutance constant
            const scalar s_;

            //- Nocturnal stomatal condutance (mol / m2 s)
            const dimensionedScalar gsn_;

            //- Maximum carboxylation capacity at 25 deg.
            const dimensionedScalar Vcmax25_;

            //- Michaelis constant for CO2 fixation
            const dimensionedScalar Kc25_;

            //- Michaelis constant for O2 fixation
            const dimensionedScalar Ko25_;

            //- Assimilation constant
            const scalar gammac_;

            //- Assimilation constant
            const scalar gammao_;

            //- Assimilation constant (1/s)
            const dimensionedScalar kco_;

            //- Assimilation constant (1/s)
            const dimensionedScalar ko_;

            //- Fraction of short-wave inside PAR band
            const scalar rPAR_;

            //- Apparent Quantum Yield 
            const scalar gammaPAR_;

            //- Number of sides stomata is present
            const scalar nEvapSides_;

            //- Absolute gravitational acceleration (m/s2)
            const dimensionedScalar gabs_;

            //- Soil-root conductance constant (-)
            const dimensionedScalar alpha_;

            


    // Model time-dependent properties

            //- Leaf water potential [Pa]
            dimensionedScalar psi_L_;

            //- Root water potential [Pa]
            dimensionedScalar psi_R_;

            //- 24h average Leaf water potential [Pa]
            scalarList psi_L24_;

            //- Marginal water use efficiency
            dimensionedScalar lambda_;

            //- Photosynthetically active radiation
            dimensionedScalar Qp_;

            //- Internal clock
            scalar internalTime;

            //- Varying vegetation properties
            IOdictionary varyingVegetationProperties_;
            


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

            //- Short-wave radiation absorbed by leaves (W/m2)
            volScalarField qrsw_leaf_;

            //- Long-wave radiation absorbed by leaves (W/m2)
            volScalarField qrlw_leaf_;

            //- Short + long wave radiation absorbed by leaves (W/m2)
            volScalarField qrad_leaf_;

            //- Stomatal condutance to CO2 (mol / m2 s)      
            volScalarField gs_;

            //- Effective conductance to CO2 (mol / m2 s)
            volScalarField gc_eff_;

            //- Effectuve conductance to vapour (mol / m2 s)
            volScalarField gv_eff_;

            //- Assimilation rate (mol / m2 s)
            volScalarField An_;

            //- Intercellular CO2 concentration (mol / m2 s)
            volScalarField ci_;            

            //- Leaf transpiration rate (kg / m2 s)
            volScalarField gv_leaf_;

            //- Latent heat flux (W / m2)
            volScalarField qlat_leaf_;

            //- Sensible heat flux (W / m2)
            volScalarField qsen_leaf_;

            //- Energy source term [kg/m s3]
            volScalarField Sh_;

            //- Momentum source term [kg/m2 s2]
            volVectorField Su_;

            //- Absolute mositure source term [kg/m3 s]
            volScalarField Sw_;                 
            
            //- CO2 source term [mol / mol s]
            volScalarField Sc_; 

            //- Soil moisture source term [kg/m3 s]
            volScalarField Sws_;

            
    // Protected Member Functions

        //- Solve the quadratic equation
        virtual scalar solve_quadratic(const scalar& qb, const scalar& qc);

        //- Calculate the intercellular CO2 concentration
        virtual scalar calc_ci(const scalar& a1, const scalar& a2, const scalar& ccp, const scalar& gc_eff, const scalar& c, const scalar& Rd);

        //- Calculate stomatal condutance (mol / m2 s)
        virtual scalar calc_gs(const scalar& a1, const scalar& a2, const scalar& VPD, const scalar& c);

        //- Calculate effective CO2 condutance (mol / m2 s)
        virtual scalar calc_gc_eff(const scalar& gs, const scalar& ga);

        //- Calculate effective vapour conductance (mol / m2 s)
        virtual scalar calc_gv_eff(const scalar& gc_eff);

        //- Calculate Assimilation rate (mol / m2 s)
        virtual scalar calc_An(const scalar& gc_eff, const scalar& ci, const scalar& c);

        //- Calculate minimum assimilation rate, stomatal resistances and intercellular CO2 concentration
        virtual void solve_assimilation(const scalar& Tl, const scalar& VPD, const scalar& c, const label& cellI);

        //- Assert the units are as required
        template<class tmpClass>
        void assertDimensions(const tmpClass& sourceVar, const dimensionSet& targetUnit);
        
        //- Write time-varying vegetation properties
        virtual void writeVegetationProperties();

        //- Return saturated water vapour pressure [Pa]
        virtual scalar calc_pvsat(const scalar& T);

        //- Return vapour pressure [Pa]
        virtual scalar calc_pv(const scalar& p, const scalar& w);

        //- Return vapour pressure deficit [Pa/Pa]
        virtual scalar calc_VPD(const scalar& T, const scalar& p, const scalar& w);

        //- Calculate aerodynamic condutance (mol/m2 s)
        virtual void calc_ga(const volVectorField& U, const volScalarField& rho);

        //- Calculate soil-root system conductance (s/m)
        virtual void calc_gsr(const volScalarField& Kl);

        //- Calculate soil xylem condutance (1/s)
        virtual scalar calc_gx(const scalar& psi_L);

        //- Calculate marginal water use efficiency
        virtual void calc_marginalWUE();

        //- calculate radiation
        virtual void calc_radiation();
        

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
            const volScalarField& c,    //- Fluid CO2 concentration [mol/mol]
			//const volScalarField& Ts,   //- Soil temperature [K]
            //const volScalarField& ws,   //- Soil moisture [kg/m3]
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
    
        //- Solve leaf energy balance
        //virtual void solve(const volVectorField&U, const rhoThermo& thermo, const volScalarField& w, const volScalarField& ws, const volScalarField& pc, const volScalarField& Kl);
        virtual void solve(const volVectorField& U, const rhoThermo& thermo, const volScalarField& w, const volScalarField& c);

        //- Solve Soil-Plant-Atmospheric Contiuum
        /////// virtual void solve_SPAC(const volScalarField& ws, const volScalarField& pc, const volScalarField& Kl);


    // Source terms

        //- Return vegetation energy source
        virtual tmp<volScalarField> Sh();

        //- Return vegetation specific humidity source
        virtual tmp<volScalarField> Sw();

        //- Return vegetation CO2 source
        virtual tmp<volScalarField> Sc();

        //- Return vegetation momentum source
        virtual tmp<fvVectorMatrix> Su(volScalarField& rho, volVectorField& U);

        //- Return soil moisture source term
        virtual tmp<volScalarField> Sws(volScalarField& Kl, volScalarField& pc); 

        //- Read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif


// ************************************************************************* //
