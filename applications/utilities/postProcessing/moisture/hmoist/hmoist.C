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
#include "emptyFvPatch.H"

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

    IOobject qheader
    (
        "q",
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

    dimensionedScalar Cpv
    (
         transportProperties.lookup("Cpv")
    );

    dimensionedScalar lambda
    (
         transportProperties.lookup("lambda")
    );

    // Fluid density
    dimensionedScalar rho
    (
        transportProperties.lookup("rho")
    );


    dimensionedScalar TRef
    (
         transportProperties.lookup("TRef")
    );

    dimensionedScalar qRef
    (
         transportProperties.lookup("qRef")
    );

    if (qheader.headerOk() && Theader.headerOk())
    {
        Info<< "    Reading q" << endl;
        volScalarField q(qheader, mesh);

        Info<< "    Reading T" << endl;
        volScalarField T(Theader, mesh);

        // specific enthalpy of dry air hda - ASHRAE 1.8
        volScalarField hda("hda", rho*Cpa*T-rho*Cpa*TRef);

        // specific enthalpy of dry vapor
        volScalarField hdv("hdv", rho*q*(lambda + Cpv*T) - rho*qRef*(lambda + Cpv*TRef));


        // specific enthalpy for moist air hmoist - ASHRAE 1.8
        volScalarField hmoist("hmoist", hda + hdv);

        if (writeResults)
        {
            hda.write();
            hdv.write();
            hmoist.write();
        }
        else
        {
            Info<< "        Min hda    : " << min(hda).value() << " [J/m3]"
                << "\n        Max hda    : "<< max(hda).value() << " [J/m3]" << endl;
            Info<< "        Min hdv    : " << min(hdv).value() << " [J/m3]"
                << "\n        Max hdv    : "<< max(hdv).value() << " [J/m3]" << endl;
            Info<< "        Min hmoist : " << min(hmoist).value() << " [J/m3]"
                << "\n        Max hmoist : "<< max(hmoist).value() << " [J/m3]" << endl;
        }


        // print results
       //
      //   surfaceScalarField hdaBoundary =
      //     fvc::interpolate(hda);
       //
      //   const surfaceScalarField::GeometricBoundaryField& patchhda =
      //       hdaBoundary.boundaryField();
      //   // const surfaceScalarField::GeometricBoundaryField& patchCondHeatFlux =
      //   //     condHeatFluxNormal.boundaryField();
      //   // const surfaceScalarField::GeometricBoundaryField& patchTurbHeatFlux =
      //   //     turbHeatFluxNormal.boundaryField();
      //   // const surfaceScalarField::GeometricBoundaryField& patchTotHeatFlux =
      //   //     totHeatFluxNormal.boundaryField();
       //
      //   Info<< "\nHeat at the boundaries " << endl;
      //   forAll(patchhda, patchi)
      //   {
      //       if ( (!isA<emptyFvPatch>(mesh.boundary()[patchi])) &&
      //            (mesh.boundary()[patchi].size() > 0) )
      //       {
      //           Info<< "    "
      //               << mesh.boundary()[patchi].name()
      //               << "\n        Total area [m2]                   : "
      //               << gSum(mesh.magSf().boundaryField()[patchi])
      //               << "\n        Integral Energy [J] : "
      //               << gSum
      //                  (
      //                      mesh.magSf().boundaryField()[patchi]
      //                     *patchhda[patchi]
      //                  )
      //               << nl << endl;
      //       }
      //  }
      //  Info << endl;

    }
    else
    {
        Info<< "    No q or No T" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
