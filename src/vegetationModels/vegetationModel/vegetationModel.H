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
    // Model coefficients
    dimensionedScalar C_; // proportionality factor of ra
    dimensionedScalar Cdf_; // leaf drag coefficient
    dimensionedScalar l_; // characteristic length of leaf
    dimensionedScalar rs_; // stomatal resistance
    dimensionedScalar Tl_; // temperature of leaf
    dimensionedScalar UMin_; //- Lower limit of U
    dimensionedScalar ql_; //- saturated humidity at leaf
    dimensionedScalar lambda_; // latent heat of vaporization
    // -----------------------------------------------
    // Model fields
    volScalarField a_; // leaf area density

public:
    // runtime type info
    TypeName("vegetationModel");

    // constructors
        vegetationModel
        (
            const volVectorField& U,
            const volScalarField& a
        );

    // destructor
        virtual ~vegetationModel()
        {}

    // Member Functions

        // -----------------------------------------------
        // return proportiality factor
        dimensionedScalar C() const
        {
            return C_;
        }

        // return leaf drag coefficient
        dimensionedScalar Cdf() const
        {
            return Cdf_;
        }

        // return characteristic length of leaf
        dimensionedScalar l() const
        {
            return l_;
        }

        // return temperature of leaf
        dimensionedScalar Tl() const
        {
            return Tl_;
        }

        // return leaf area density
        tmp<volScalarField> a() const
        {
            return a_;
        }

        // return vegetation properties dictionary
        dictionary vegetationProperties() const
        {
            return vegetationProperties_;
        }


        // -----------------------------------------------

        // return aerodynamic resistance
        virtual tmp<volScalarField> ra(volVectorField& U) const;
        // return sensible heat
        virtual tmp<fvScalarMatrix> Qs(volVectorField& U, volScalarField& T) const;

        // return tranpiration rate
        virtual tmp<fvScalarMatrix> E(volVectorField& U, volScalarField& T, volScalarField& q) const;

        // return latent heat
        virtual tmp<fvScalarMatrix> Ql(volVectorField& U, volScalarField& T, volScalarField& q) const;


        // -----------------------------------------------
        // return vegetation momentum source
        virtual tmp<fvVectorMatrix> Su(volVectorField& U) const;

        // return vegetation energy source
        virtual tmp<fvScalarMatrix> Sh(volVectorField& U, volScalarField& T) const;

        // return vegetation specific humidity source
        virtual tmp<fvScalarMatrix> Sq(volVectorField& U, volScalarField& T, volScalarField& q) const;


        // -----------------------------------------------
        // read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif

// ************************************************************************* //
