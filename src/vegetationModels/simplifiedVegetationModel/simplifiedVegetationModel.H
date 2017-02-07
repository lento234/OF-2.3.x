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
    Foam::simplifiedVegetationModel

Description
    Simplified Vegetation model

Author
    Lento Manickathan. manickathan@arch.ethz.ch

SourceFiles
    vegetationModel.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedVegetationModel_H
#define simplifiedVegetationModel_H

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                         Class vegetationModel Declaration
\*---------------------------------------------------------------------------*/

class simplifiedVegetationModel
:
    public IOdictionary

{

private:
    // disallow copy construct
    simplifiedVegetationModel(const simplifiedVegetationModel&);

    // disallow default bitwise assignment
    void operator=(const simplifiedVegetationModel&);

protected:
    // -----------------------------------------------
    dictionary vegetationProperties_;
    // -----------------------------------------------
    const Time& runTime_;       // time
    const fvMesh& mesh_;        // mesh

    // -----------------------------------------------
    // Model coefficients
    dimensionedScalar cpa_;     // specific heat of air at constant pressure
    dimensionedScalar E_;    // transpiration rate
    dimensionedScalar Qs_;    // ensible heat flux
    dimensionedScalar rhoa_;    // denisty of air

    // -----------------------------------------------
    // Model fields
    volScalarField Cf_;     // vegetation drag coefficient
    volScalarField Sh_;     // energy source term
    volScalarField Sq_;     // humidity source term
    volVectorField Su_;     // momentum source term



public:
    // runtime type info
    TypeName("simplifiedVegetationModel");

    // constructors
        simplifiedVegetationModel
        (
            const volVectorField& U,
            const volScalarField& T,
            const volScalarField& q
        );

    // destructor
        virtual ~simplifiedVegetationModel()
        {}

    // Member Functions

        // -----------------------------------------------

        // return vegetation properties dictionary
        dictionary vegetationProperties() const
        {
            return vegetationProperties_;
        }

        // -----------------------------------------------


        // solve all
        virtual void solve(volVectorField&U, volScalarField& T, volScalarField& q);

        // -----------------------------------------------
        // source terms

        // return vegetation energy source
        virtual tmp<volScalarField> Sh();

        // return vegetation momentum source
        virtual tmp<fvVectorMatrix> Su(volVectorField& U);

        // return vegetation specific humidity source
        virtual tmp<volScalarField> Sq();

        // -----------------------------------------------

        // read vegetationProperties dictionary
        virtual bool read();
};


} // end namepsace Foam

#endif


// ************************************************************************* //
