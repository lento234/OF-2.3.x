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
    Foam::vegetationModel

Description
    Vegetation model

Author
    Lento Manickathan. manickathan@arch.ethz.ch

SourceFiles
    vegetationModel.C

\*---------------------------------------------------------------------------*/

#ifndef vegetationModel_H
#define vegetationModel_H

#include <typeinfo>
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                         Class vegetationModel Declaration
\*---------------------------------------------------------------------------*/

class vegetationModel
:
    public IOdictionary

{

private:
    // disallow copy construct
    vegetationModel(const vegetationModel&);

    // disallow default bitwise assignment
    void operator=(const vegetationModel&);

protected:
    // -----------------------------------------------
    dictionary vegetationProperties_;
    // -----------------------------------------------
    const Time& runTime_;       // time
    const fvMesh& mesh_;        // mesh
    // -----------------------------------------------
    // Model coefficients
    dimensionedScalar cpa_;     // specific heat of air at constant pressure
    dimensionedScalar C_;       // proportionality factor of ra
    dimensionedScalar Cdf_;     // leaf drag coefficient
    dimensionedScalar H_;       // Height of vegetation
    dimensionedScalar kc_;      // extinction coefficient of radiation
    dimensionedScalar l_;       // characteristic length of leaf
    dimensionedScalar Rg0_;     // global radiation at top of canopy
    dimensionedScalar rhoa_;    // denisty of air
    dimensionedScalar rs_;      // stomatal resistance
    dimensionedScalar TlMin_;   // lower limit of T
    dimensionedScalar UMin_;    // lower limit of U
    dimensionedScalar lambda_;  // latent heat of vaporization
    // -----------------------------------------------
    // Model fields
    volScalarField LAD_;    // leaf area density
    volScalarField LAI_;    // leaf area index
    volScalarField E_;      // transpiration rate
    volScalarField qsat_;   // saturated specific humidity
    volScalarField Ql_;     // latent heat flux
    volScalarField ra_;     // aerodynamic resistance
    volScalarField rhosat_; // saturated density field
    volVectorField Rg_;     // global radiation
    volScalarField Rn_;     // global radiation density in volume
    volScalarField Sq_;     // humidity source term
    volVectorField Su_;     // momentum source term
    volScalarField Tl_;     // temperature of leaf


public:
    // runtime type info
    TypeName("vegetationModel");

    // constructors
        vegetationModel
        (
            const volVectorField& U,
            const volScalarField& LAD,
            const volScalarField& LAI,
            const volScalarField& T
        );

    // destructor
        virtual ~vegetationModel()
        {}

    // Member Functions

        // -----------------------------------------------

        // return proportiality factor
        // dimensionedScalar C() const
        // {
        //     return C_;
        // }

        // return leaf drag coefficient
        // dimensionedScalar Cdf() const
        // {
        //     return Cdf_;
        // }

        // return characteristic length of leaf
        // dimensionedScalar l() const
        // {
        //     return l_;
        // }

        // return leaf area density
        // tmp<volScalarField> LAD() const
        // {
        //     return LAD_;
        // }

        // return radiation density inside vegetation
        tmp<volScalarField> Rn() const
        {
            return Rn_;
        }

        // return leaf area density
        // tmp<volScalarField> LAD() const
        // {
        //     return LAD_;
        // }

        // dimensionedScalar rs() const
        // {
        //     return rs_;
        // }

        // return vegetation properties dictionary
        dictionary vegetationProperties() const
        {
            return vegetationProperties_;
        }

        // -----------------------------------------------

        // return transpiration rate
        // virtual tmp<volScalarField> E(volVectorField& U, volScalarField& T, volScalarField& q) const;

        // return latent heat
        // virtual tmp<volScalarField> Ql(volVectorField& U, volScalarField& T, volScalarField& q) const;

        // return sensible heat
        // virtual tmp<volScalarField> Qs(volVectorField& U, volScalarField& T) const;



        // return saturated density of water vapour
        virtual double calc_rhosat(double& T);
        // virtual tmp<volScalarField> rhosat(volScalarField& T) const;

        // return empirical global radiation
        //virtual tmp<volScalarField> Rg(volVectorField& U) const;
        // virtual tmp<volVectorField> Rg();

        // return net radiation density in volume
        // virtual tmp<volScalarField> Rn();

        // return saturated specific humidity of water vapour
        // virtual double calc_qsat(double& rhosat);
        // virtual tmp<volScalarField> qsat(volScalarField& T) const;


        // -----------------------------------------------
        // source terms

        // return vegetation energy source
        // virtual tmp<fvScalarMatrix> Sh(volVectorField& U, volScalarField& T) const;

        // return vegetation momentum source
        // virtual tmp<fvVectorMatrix> Su(volVectorField& U) const;
        virtual tmp<volVectorField> Su(volVectorField& U);

        // return vegetation specific humidity source
        virtual tmp<volScalarField> Sq();
        // virtual tmp<fvScalarMatrix> Sq(volVectorField& U, volScalarField& T, volScalarField& q) const;

        // -----------------------------------------------

        // calc radiation
        virtual void radiation();

        // calc aerodynamic, stomatal resistances
        virtual void resistance(volVectorField& U);

        // solve all
        virtual void solve(volVectorField&U, volScalarField& T, volScalarField& q);

        // testing functions
        // virtual void testing(volVectorField& U, volScalarField& T, volScalarField& q);


        // read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif

// ************************************************************************* //
