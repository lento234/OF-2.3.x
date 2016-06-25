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

    IOobject Wheader
    (
        "W",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject Theader
    (
        "T",
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

    dimensionedScalar p
    (
         transportProperties.lookup("p")
    );


    if (Wheader.headerOk() && Theader.headerOk())
    {
        Info<< "    Reading W" << endl;
        volScalarField W(Wheader, mesh);

        Info<< "    Reading T" << endl;
        volScalarField T(Theader, mesh);


        // saturated vapor pressure pws - ASHRAE 1.2
        dimensionedScalar c0("c0", dimPressure, 1);
        dimensionedScalar c1("c1", dimTemperature, -5.8002206e3);
        dimensionedScalar c2("c2", dimless, 1.3914993);
        dimensionedScalar c3("c3", dimless/dimTemperature, -4.8640239e-2);
        dimensionedScalar c4("c4", dimless/dimTemperature/dimTemperature, 4.1764768e-5);
        dimensionedScalar c5("c5", dimless/dimTemperature/dimTemperature/dimTemperature, -1.4452093e-8);
        dimensionedScalar c6("c6", dimless, 6.5459673);
        dimensionedScalar c7("c7", dimless/dimTemperature, 1);

        volScalarField pws("pws", c0*exp(  c1/T
                                         + c2
                                         + c3*T
                                         + c4*pow(T,2)
                                         + c5*pow(T,3)
                                         + c6*log(c7*T) ));

        // vapor pressure pw - ASHRAE 1.8
        volScalarField pw("pw", p*W/(0.621945+W));

        Info<< "    Calculating RH" << endl;
        volScalarField RH("RH", (pw/pws)*100);

        if (writeResults)
        {
            RH.write();
        }
        else
        {
            Info<< "        Min W : " << min(RH).value() << " [kg_w/kg_da]"
                << "\n        Max W : "<< max(RH).value() << " [kg_w/kg_da]" << endl;
        }
    }
    else
    {
        Info<< "    No W or No T" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
