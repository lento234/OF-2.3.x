/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    YvToq

Description
    Calculates and writes the specific humdity q =
    (mass of water vapor) / (mass of dry air) from the mass fraction of
    water vapour Yv = (mass of water vapour) / (mass of moist air)

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Yvheader
    (
        "Yv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Yvheader.headerOk())
    {
        Info<< "    Reading Yv" << endl;
        volScalarField Yv(Yvheader, mesh);

        Info<< "    Calculating q" << endl;
        volScalarField q
        (
            IOobject
            (
                "q",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            Yv/(1.0-Yv)
        );

        if (writeResults)
        {
            q.write();
        }
        else
        {
            Info<< "        Min q : " << gMin(q) << " [-]"
                << "\n        Max q : "<< gMax(q) << " [-]" << endl;
        }
    }
    else
    {
        Info<< "    No q" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
