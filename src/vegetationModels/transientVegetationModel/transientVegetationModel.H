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
    Foam::transientVegetationModel

Description
    Vegetation model

Author
    Lento Manickathan. manickathan@arch.ethz.ch

SourceFiles
    transientVegetationModel.C

\*---------------------------------------------------------------------------*/

#ifndef transientVegetationModel_H
#define transientVegetationModel_H

// #include <typeinfo>
// #include "RASModel.H"
// #include "incompressible/turbulenceModel/turbulenceModel.H"
// #include "volFields.H"
// #include "surfaceFields.H"
#include "fvm.H"
// #include "fvc.H"
// #include "fvMatrices.H"
// #include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
// #include "Switch.H"
#include "bound.H"
// #include "autoPtr.H"
// #include "runTimeSelectionTables.H"
#include "TimeDataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                         Class transientVegetationModel Declaration
\*---------------------------------------------------------------------------*/

class transientVegetationModel
:
    public IOdictionary

{

private:
    // disallow copy construct
    transientVegetationModel(const transientVegetationModel&);

    // disallow default bitwise assignment
    void operator=(const transientVegetationModel&);

protected:
    // -----------------------------------------------
    dictionary vegetationProperties_;
    // -----------------------------------------------
    const Time& runTime_;       // time
    const fvMesh& mesh_;        // mesh
    // -----------------------------------------------
    // Model coefficients
    dimensionedScalar a1_;     // stomatal resistance specific model parameters
    dimensionedScalar a2_;     // stomatal resistance specific model parameters
    dimensionedScalar a3_;     // stomatal resistance specific model parameters
    dimensionedScalar cpa_;     // specific heat of air at constant pressure
    dimensionedScalar C_;       // proportionality factor of ra
    dimensionedScalar D0_;      // vapour pressure corresponding to minimal resistance
    dimensionedScalar H_;       // Height of vegetation
    dimensionedScalar kc_;      // extinction coefficient of radiation
    dimensionedScalar l_;       // characteristic length of leaf
    TimeDataEntry<scalar> Rg0_;     // global radiation at top of canopy
    TimeDataEntry<scalar> Rl0_;     // long-wave radiation at top of canopy
    dimensionedScalar rhoa_;    // denisty of air
    dimensionedScalar rsMin_;   // minimum stomatal resistance
    dimensionedScalar TlMin_;   // lower limit of T
    dimensionedScalar UMin_;    // lower limit of U
    dimensionedScalar lambda_;  // latent heat of vaporization
    // -----------------------------------------------
    // Model fields
    volScalarField Cf_;     // vegetation drag coefficient
    volScalarField E_;      // transpiration rate
    volScalarField pv_;     // water vapor pressure
    volScalarField pvsat_;  // saturated water vapor pressure
    volScalarField LAD_;    // leaf area density
    volScalarField LAI_;    // leaf area index
    volScalarField Yvsat_;   // saturated specific humidity
    volScalarField Ql_;     // latent heat flux
    volScalarField Qs_;     // sensible heat flux
    volScalarField ra_;     // aerodynamic resistance
    volScalarField rs_;     // stomatal resistance
    volScalarField rhosat_; // saturated density field
    volVectorField Rg_;     // global radiation
    volScalarField Rn_;     // global radiation density in volume
    volScalarField Sh_;     // energy source term
    volScalarField SYv_;     // humidity source term
    volVectorField Su_;     // momentum source term
    volScalarField Tl_;     // temperature of leaf
    volScalarField VPD_;    // vapor pressure deficit

    double cpa, C, H, kc, l, lambda, Rg0, rhoa, Rl0, rsMin;

public:
    // runtime type info
    TypeName("transientVegetationModel");

    // constructors
        transientVegetationModel
        (
            const volVectorField& U,
            const volScalarField& T,
            const volScalarField& Yv
        );

    // destructor
        virtual ~transientVegetationModel()
        {}

    // Member Functions

        // -----------------------------------------------

        // return vegetation properties dictionary
        dictionary vegetationProperties() const
        {
            return vegetationProperties_;
        }

        // -----------------------------------------------

        // return saturated density of water vapour
        virtual double calc_rhosat(double& T);

        // return saturated water vapour pressure
        virtual double calc_pvsat(double& T);

        // calc radiation
        virtual void radiation();

        // calc aerodynamic, stomatal resistances
        virtual void resistance(volScalarField& magU, volScalarField& T, volScalarField& Yv);

        // solve all
        virtual void solve(volVectorField&U, volScalarField& T, volScalarField& Yv);

        // -----------------------------------------------
        // source terms

        // return vegetation energy source
        virtual tmp<volScalarField> Sh();

        // return vegetation momentum source
        virtual tmp<volVectorField> Su(volVectorField& U);

        // return vegetation specific humidity source
        virtual tmp<volScalarField> SYv();

        // -----------------------------------------------

        // read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif

// return proportiality factor
// dimensionedScalar C() const
// {
//     return C_;
// }

// ************************************************************************* //
