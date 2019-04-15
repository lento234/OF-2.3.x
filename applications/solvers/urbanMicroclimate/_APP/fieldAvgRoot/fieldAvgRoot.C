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
    fieldAvg

Description
    Calculates the average value of specific field at the vegetation zone

Author
    Lento Manickathan <lento.manickathan@gmail.com>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
void printIntegrate
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const IOobject& RADFieldHeader,
    bool& done
)
{
    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        Info<< "    Reading " << fieldHeader.headerClassName() << " "
            << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);
        volScalarField RAD(RADFieldHeader, mesh);

        Info<< "    Integral of " << fieldHeader.name()
            << " over full fluid volume = "
            << gSum(mesh.V()*field.internalField()*RAD.internalField())/gSum(mesh.V()*RAD.internalField()) << " "
            << field.dimensions()
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

    // Give fluid volume
    Info<< "    Volume of fluid = "
        << gSum(mesh.V()) << " [m3]" << endl;


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

        IOobject RADFieldHeader
        (
            "RAD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        // Check field exists
        if (fieldHeader.headerOk() && RADFieldHeader.headerOk())
        {
            mesh.readUpdate();
            // Read field and calc integral
            bool done = false;

            printIntegrate<volScalarField>(mesh,fieldHeader,RADFieldHeader,done);
            printIntegrate<volVectorField>(mesh,fieldHeader,RADFieldHeader,done);
            printIntegrate<volSphericalTensorField>(mesh,fieldHeader,RADFieldHeader,done);
            printIntegrate<volSymmTensorField>(mesh,fieldHeader,RADFieldHeader,done);
            printIntegrate<volTensorField>(mesh,fieldHeader,RADFieldHeader,done);

            if (!done)
            {
                FatalError
                    << "Only possible to integrate volFields."
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
