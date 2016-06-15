/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    cellIntegrate

Description
    Calculates the integral of the specified field over the specified patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
void printIntegrate
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    bool& done
)
{
    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        Info<< "    Reading " << fieldHeader.headerClassName() << " "
            << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);

        Info<< "    Integral of " << fieldHeader.name()
            << " over full fluid volume = "
            << gSum(mesh.V()*field.internalField())
            << nl;

        done = true;
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::validArgs.append("fieldName");
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const word fieldName = args[1];

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field exists
        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            // Give fluid volume
            Info<< "    Volume of fluid = "
                << gSum(mesh.V()) << endl;

            // Read field and calc integral
            bool done = false;

            printIntegrate<volScalarField>
            (
                mesh,
                fieldHeader,
                done
            );
            printIntegrate<volVectorField>
            (
                mesh,
                fieldHeader,
                done
            );

            printIntegrate<volSphericalTensorField>
            (
               mesh,
               fieldHeader,
               done
            );
            printIntegrate<volSymmTensorField>
            (
               mesh,
               fieldHeader,
               done
            );
            printIntegrate<volTensorField>
            (
               mesh,
               fieldHeader,
               done
            );

            if (!done)
            {
                FatalError
                    << "Only possible to integrate "
                    << "volFields and surfaceFields."
                    << " Field " << fieldName << " is of type "
                    << fieldHeader.headerClassName()
                    << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
