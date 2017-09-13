/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    radiationCoupledBase

Description
    Common functions to emissivity. It gets supplied from lookup into a
    dictionary or calculated by the solidThermo:

    - 'lookup' : lookup volScalarField with name
    - 'solidThermo' : use solidThermo emissivity()

SourceFiles
    radiationCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef solarRadiationCoupledBase_H
#define solarRadiationCoupledBase_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solarRadiationCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class solarRadiationCoupledBase
{
public:

    //- Type of supplied albedo
    enum albedoMethodType
    {
        SOLIDRADIATION,
        LOOKUP
    };


private:

    // Private data

        static const NamedEnum<albedoMethodType, 2>
            albedoMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;


protected:

    // Protected data

        //- How to get albedo
        const albedoMethodType method_;

         //- albedo
         // Cached locally when is read from dictionary (lookup mode)
        scalarField albedo_;


public:

    // Constructors

        //- Construct from patch, albedo mode and albedo
        solarRadiationCoupledBase
        (
            const fvPatch& patch,
            const word& calculationMethod,
            const scalarField& albedo
        );

        //- Construct from patch and dictionary
        solarRadiationCoupledBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );


    // Member functions

        //- Method to obtain albedo
        word albedoMethod() const
        {
            return albedoMethodTypeNames_[method_];
        }


        //- Calculate corresponding albedo field
        scalarField albedo() const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
