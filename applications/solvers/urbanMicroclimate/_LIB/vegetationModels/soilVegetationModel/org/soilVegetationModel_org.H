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
#include "fvm.H"
#include "IOdictionary.H"
#include "bound.H"
#include "scalarListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"

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
    
    
            const Time& runTime_;       // time
            const fvMesh& mesh_;        // mesh

        // Model coefficients
    
            dimensionedScalar a1_;          // stomatal resistance specific model parameters
            dimensionedScalar a2_;          // stomatal resistance specific model parameters
            dimensionedScalar a3_;          // stomatal resistance specific model parameters
            dimensionedScalar cpa_;         // specific heat of air at constant pressure
            dimensionedScalar C_;           // proportionality factor of ra
            dimensionedScalar D0_;          // vapour pressure corresponding to minimal resistance
            dimensionedScalar nEvapSides_;  // number of sides, leaf evaporates from
            dimensionedScalar H_;           // Height of vegetation
            dimensionedScalar kc_;          // extinction coefficient of radiation
            dimensionedScalar l_;           // characteristic length of leaf
            dimensionedScalar Rg0_;         // global radiation at top of canopy
            dimensionedScalar Rl0_;         // long-wave radiation at top of canopy
            dimensionedScalar rhoa_;        // denisty of air
            dimensionedScalar rsMin_;       // minimum stomatal resistance
            dimensionedScalar TlMin_;       // lower limit of T
            dimensionedScalar UMin_;        // lower limit of U
            dimensionedScalar lambda_;      // latent heat of vaporization

            scalarListIOList divqrsw;
    

        // Fields

            volScalarField Cf_;             // vegetation drag coefficient
            volScalarField E_;              // transpiration rate
            //volScalarField pv_;             // water vapor pressure
            //volScalarField pvsat_;          // saturated water vapor pressure
            volScalarField LAD_;            // leaf area density
            volScalarField LAI_;            // leaf area index
            //volScalarField wsat_;           // saturated specific humidity
            volScalarField Qlat_;           // latent heat flux
            volScalarField Qsen_;           // sensible heat flux
            volScalarField ra_;             // aerodynamic resistance
            volScalarField rs_;             // stomatal resistance
            //volScalarField rhosat_;         // saturated density field
            volVectorField Rg_;             // global radiation
            volScalarField Rn_;             // global radiation density in volume
            volScalarField Tl_;             // temperature of leaf
            
            volScalarField Sh_;             // energy source term
            volScalarField Sw_;             // humidity source term
            volVectorField Su_;             // momentum source term
            // volScalarField Sws_;            // soil moisture content source term
            
            //volScalarField VPD_;            // vapor pressure deficit

            

public:
    
    //-Runtime type information
    TypeName("soilVegetationModel");

    // Constructors
    
        //- Construct from components
        soilVegetationModel
        (
            const volVectorField& U,
            const volScalarField& T,
            const volScalarField& w
			//volScalarField& Tl_
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

        
        //- Return saturated density of water vapour
        virtual double calc_rhosat(double& T);

        //- Return saturated water vapour pressure
        virtual double calc_pvsat(double& T);

        //- Return vapour pressure
        virtual double calc_pv(double& p, double& w);

        //- calculate radiation
        virtual void radiation();

        //- calculate aerodynamic, stomatal resistances
        virtual void resistance(volScalarField& magU, volScalarField& T, volScalarField& w, volScalarField& Tl);

        //- Solve leaf energy balance
        virtual void solve(volVectorField&U, volScalarField& T, volScalarField& w);//, volScalarField& Tl_);

        // source terms

        //- Return vegetation energy source
        virtual tmp<volScalarField> Sh();

        //- Return vegetation specific humidity source
        virtual tmp<volScalarField> Sw();

        //- Return vegetation momentum source
        virtual tmp<fvVectorMatrix> Su(volScalarField& rho, volVectorField& U);

        // //- Return soil moisture source term
        // virtual tmp<volScalarField> Sws(); 

        //- Return net mass transpiration rate
        //dimensionedScalar mtrans();

        //- Read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif


// ************************************************************************* //
