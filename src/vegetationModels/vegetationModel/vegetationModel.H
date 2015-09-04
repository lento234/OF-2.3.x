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

SourceFiles
    vegetationModelI.H
    vegetationModel.C
    vegetationModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef vegetationModel_H
#define vegetationModel_H

#include "RASModel.H"
// #include "incompressible/turbulenceModel/turbulenceModel.H"
// #include "volFields.H"
// #include "surfaceFields.H"
// #include "fvm.H"
// #include "fvc.H"
// #include "fvMatrices.H"
// #include "incompressible/transportModel/transportModel.H"
// #include "IOdictionary.H"
// #include "Switch.H"
// #include "bound.H"
// #include "autoPtr.H"
// #include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                         Class vegetationModel Declaration
\*---------------------------------------------------------------------------*/

class vegetationModel
{

    protected:
        // Model coefficients and Fields
        volScalarField a_;
        dimensionedScalar Cdf_;

    public:
        // Constructor
        vegetationModel
        (
            const volScalarField& a,
            const dimensionedScalar& Cdf
        );

        // Destructor
        virtual ~vegetationModel()
        {}

        // Member Functions
        tmp<volScalarField> a() const
        {
            return a_;
        }

        dimensionedScalar Cdf() const
        {
            return Cdf_;
        }

        // Momentum source
        virtual tmp<fvVectorMatrix> Su(volVectorField& U) const;


};


}

#endif

// ************************************************************************* //
