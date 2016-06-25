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
    RHtoW

Description
    Convert Relative Humidity RH (%) to Humidity ratio W (kg_w/kg_da)
    Humidity ratio W is ratio of mass of water vapor to mass of dry air

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Theader
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject Wheader
    (
        "W",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Cpa
    (
         transportProperties.lookup("Cpa")
    );

    dimensionedScalar Cpw
    (
         transportProperties.lookup("Cpw")
    );

    dimensionedScalar lambda
    (
         transportProperties.lookup("lambda")
    );

    dimensionedScalar TRef
    (
         transportProperties.lookup("TRef")
    );

    dimensionedScalar WRef
    (
         transportProperties.lookup("WRef")
    );

    if (Wheader.headerOk() && Theader.headerOk())
    {
        Info<< "    Reading W" << endl;
        volScalarField W(Wheader, mesh);

        Info<< "    Reading T" << endl;
        volScalarField T(Theader, mesh);

        // specific enthalpy of dry air hda - ASHRAE 1.8
        volScalarField hda("hda", Cpa*(T-TRef));

        // specific enthalpy for moist air hmoist - ASHRAE 1.8
        volScalarField hmoist("hmoist", hda + (W-WRef)*(lambda + Cpw*(T-TRef)));

        if (writeResults)
        {
            hda.write();
            hmoist.write();
        }
        else
        {
            Info<< "        Min hda    : " << min(hda).value() << " [kg_w/kg_da]"
                << "\n        Max hda    : "<< max(hda).value() << " [kg_w/kg_da]" << endl;
            Info<< "        Min hmoist : " << min(hmoist).value() << " [kg_w/kg_da]"
                << "\n        Max hmoist : "<< max(hmoist).value() << " [kg_w/kg_da]" << endl;
        }
    }
    else
    {
        Info<< "    No W or No T" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
