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


    Lento Manickathan

\*---------------------------------------------------------------------------*/

#include "vegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vegetationModel, 0);
    defineRunTimeSelectionTable(vegetationModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

IOobject Foam::radiation::radiationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "vegetationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vegetationModel::vegetationModel
(
    const word& name,
    const dictionary& vegetationDict
)
:
        IOdictionary
        (
            IOobject
            (
                "vegetationProperties",
                LAD.time().constant(),
                LAD.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dict
        ),
        mesh_(LAD.mesh()),
        runTime_(LAD.time())

        //- Reference to time database
        const Time& runTime_;
    
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // end namespace Foam

