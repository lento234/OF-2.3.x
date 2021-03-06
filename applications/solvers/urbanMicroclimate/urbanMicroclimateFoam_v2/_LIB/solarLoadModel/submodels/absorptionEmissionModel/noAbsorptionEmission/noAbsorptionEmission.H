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
    Foam::radiation::noAbsorptionEmission

Description
    Dummy absorption-emission model for 'none'

SourceFiles
    noAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef solarLoadNoAbsorptionEmission_H
#define solarLoadNoAbsorptionEmission_H

#include "solarLoadAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solarLoad
{

/*---------------------------------------------------------------------------*\
                  Class noAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class noAbsorptionEmission
:
    public solarLoadAbsorptionEmissionModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~noAbsorptionEmission();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solarLoad
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
