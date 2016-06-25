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
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// saturated vapor pressure pws
volScalarField calc_pws(volScalarField& T)
{
    // saturated vapor pressure pws - ASHRAE 1.2
    dimensionedScalar c0("c0", dimPressure, 1);
    dimensionedScalar c1("c1", dimTemperature, -5.8002206e3);
    dimensionedScalar c2("c2", dimless, 1.3914993);
    dimensionedScalar c3("c3", dimless/dimTemperature, -4.8640239e-2);
    dimensionedScalar c4("c4", dimless/dimTemperature/dimTemperature, 4.1764768e-5);
    dimensionedScalar c5("c5", dimless/dimTemperature/dimTemperature/dimTemperature, -1.4452093e-8);
    dimensionedScalar c6("c6", dimless, 6.5459673);
    dimensionedScalar c7("c7", dimless/dimTemperature, 1);

    return c0*exp(  c1/T + c2 + c3*T + c4*pow(T,2) + c5*pow(T,3) + c6*log(c7*T) );

}

// vapor pressure pw
volScalarField calc_pw(volScalarField& RH, volScalarField& T)
{
     // saturated vapor pressure pws
     volScalarField pws("pws", calc_pws(T));

     // vapor pressure pw - ASHRAE 1.8
     return (RH/100.0)*pws;
}

// humidity ratio W
volScalarField calc_W(volScalarField& RH, volScalarField& T, dimensionedScalar& p)
{
      // vapor pressure pw
      volScalarField pw("pw", calc_pw(RH, T));

     // humidity ratio
     return 0.621945*pw/(p-pw);
}

// humidity ratio W from T, Twb and p
volScalarField calc_W_method2(volScalarField& T, volScalarField& Twb, dimensionedScalar& p)
{

  volScalarField pws("pws", calc_pws(Twb));

  volScalarField Ws("Ws", 0.621945 * pws/(p - pws));

  dimensionedScalar TK0("TK0", dimTemperature, -273.15);

  volScalarField t("t", T+TK0);
  volScalarField twb("twb", Twb+TK0);

  dimensionedScalar lambda("lambda", dimTemperature, 2501);

  return (((lambda - 2.326*twb)*Ws - 1.006*(t - twb))/
                (lambda + 1.86*t - 4.186*twb));
}

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

    IOobject RHheader
    (
        "RH",
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


    if (RHheader.headerOk() && Theader.headerOk())
    {
        Info<< "    Reading RH" << endl;
        volScalarField RH(RHheader, mesh);

        Info<< "    Reading T" << endl;
        volScalarField T(Theader, mesh);

        // calculate wet-bulb
        volScalarField W_target("W_target", calc_W(RH, T, p));

        // initial guess
        volScalarField Twb("Twb", T);

        // calculate W_guess
        volScalarField W_guess("W_guess", calc_W_method2(T, Twb, p));

        dimensionedScalar Teps("Teps", dimTemperature, 0.001);

        int i=0;
        while (max(mag(W_guess-W_target)/mag(W_target)).value() > 0.0001)
        {
             Info << "error" << min(mag(W_guess-W_target)) << endl;

             // calculate W_guess
             volScalarField W_guess("W_guess", calc_W_method2(T, Twb, p));

             // determine gradient
             volScalarField Twb_eps("Twb_eps", Twb-Teps);
             volScalarField W_guess_eps("W_guess2", calc_W_method2(T, Twb_eps, p));

             // gradient
             volScalarField dW_dTwb("dW_dTwb", (W_guess-W_guess_eps)/Teps);

            //  Info << "diff: " << max((W_guess-W_target)/dW_dTwb) << endl;

             // Newton-Rhapson iterative guess
             Twb -= (W_guess-W_target)/dW_dTwb;

             i++;
             if (i>1000)
                 break;
        }

        // specific enthalpy of dry air hda - ASHRAE 1.8
        // volScalarField Twb("Twb", calc_W(RH, T, p));


        if (writeResults)
        {
            Twb.write();
        }
        else
        {
            Info<< "        Min Twb    : " << min(Twb).value() << " [kg_w/kg_da]"
                << "\n        Max Twb    : "<< max(Twb).value() << " [kg_w/kg_da]" << endl;
        }
    }
    else
    {
        Info<< "    No RH or No T" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
