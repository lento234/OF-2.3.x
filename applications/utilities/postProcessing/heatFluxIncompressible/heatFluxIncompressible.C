/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    heatFluxIncompressible

Description
    Calculates and writes the all the heat flux for fluid domain and
    with (print) option, prints the fluxes at the boundaries.
    Heat fluxes:
      - (mean) convective heat flux qconv = rho*Cp*U*(T-Tref)
      - laminar conductive heat flux qcond = rho*Cp*alpha*gradT
      - turbulent convective heat flux qturb = rho*Cp*alphat*gradT
      - total heat flux qtot = qconv + qcond + qturb

    Prints the integrated flux for all (non-empty) patches.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    argList::addBoolOption
    (
      "print",
      "print heat fluxes at the boundaries"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const bool printResults = args.optionFound("print");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"
        #include "readTransportProperties.H"

         // update the turbulence fields
        turbulence->read();

        if (!(IOobject("alpha", runTime.timeName(), mesh).headerOk()))
        {
            Info<< "\nCalculating laminar heat conductivity " << endl;
            alpha = turbulence->nu()/Pr;
            alpha.correctBoundaryConditions();
        }
        else
        {
            Info<< "\nRead turbulent heat conductivity alpha" << endl;
        }

        if (!(IOobject("alphat", runTime.timeName(), mesh).headerOk()))
        {
            Info<< "\nCalculating turbulent heat conductivity " << endl;
            alphat = turbulence->nut()/Prt;
            alphat.correctBoundaryConditions();
        }
        else
        {
            Info<< "\nRead turbulent heat conductivity alphat" << endl;
        }

        if (!(IOobject("alphaEff", runTime.timeName(), mesh).headerOk()))
        {
            Info<< "\nCalculating effective heat conductivity " << endl;
            alphaEff=alpha+alphat;
            alphaEff.correctBoundaryConditions();
        }
        else
        {
            Info<< "\nRead effective heat conductivity alphaEff" << endl;
        }

        // calculate mean convective heat flux
        volVectorField convHeatFlux("convHeatFlux", rho*Cp*U*(T-TRef));
        Info<< "\nWriting convHeatFlux field" << endl;    
        convHeatFlux.write();


        // calculate gradient of temperature
        volVectorField gradT("gradT", fvc::grad(T));
        gradT.correctBoundaryConditions();

        // calculate laminar heat conductivity
        volVectorField condHeatFlux("condHeatFlux", rho*Cp*alpha*gradT);

        // calculate turbulent convective heat flux
        volVectorField turbHeatFlux("turbHeatFlux", rho*Cp*alphat*gradT);

        // calculate total heat flux
        volVectorField totHeatFlux("totHeatFlux", convHeatFlux+condHeatFlux+turbHeatFlux);

        Info<< "\nWriting condHeatFlux field" << endl;
        condHeatFlux.write();

        Info<< "\nWriting turbHeatFlux field" << endl;
        turbHeatFlux.write();

        Info<< "\nWriting totHeatFlux field" << nl << endl;
        totHeatFlux.write();

        // print results
        if (printResults)
        {

            surfaceScalarField convHeatFluxNormal =
              fvc::interpolate(convHeatFlux) & mesh.Sf()/mesh.magSf();
            surfaceScalarField condHeatFluxNormal =
              fvc::interpolate(condHeatFlux) & mesh.Sf()/mesh.magSf();
            surfaceScalarField turbHeatFluxNormal =
              fvc::interpolate(turbHeatFlux) & mesh.Sf()/mesh.magSf();
            surfaceScalarField totHeatFluxNormal =
              fvc::interpolate(totHeatFlux) & mesh.Sf()/mesh.magSf();

            const surfaceScalarField::GeometricBoundaryField& patchConvHeatFlux =
                convHeatFluxNormal.boundaryField();
            const surfaceScalarField::GeometricBoundaryField& patchCondHeatFlux =
                condHeatFluxNormal.boundaryField();
            const surfaceScalarField::GeometricBoundaryField& patchTurbHeatFlux =
                turbHeatFluxNormal.boundaryField();
            const surfaceScalarField::GeometricBoundaryField& patchTotHeatFlux =
                totHeatFluxNormal.boundaryField();

            Info<< "\nHeat fluxes at the boundaries with TRef [K]: " << TRef.value() << endl;
            forAll(patchConvHeatFlux, patchi)
            {
                if ( (!isA<emptyFvPatch>(mesh.boundary()[patchi])) &&
                     (mesh.boundary()[patchi].size() > 0) )
                {
                    Info<< "    "
                        << mesh.boundary()[patchi].name()
                        << "\n        Total area [m2]                   : "
                        << gSum(mesh.magSf().boundaryField()[patchi])
                        << "\n        Integral convective heat flux [W] : "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchConvHeatFlux[patchi]
                           )
                        << "\n        Integral conductive heat flux [W] : "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchCondHeatFlux[patchi]
                           )
                        << "\n        Integral turbulent heat flux [W]  : "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchTurbHeatFlux[patchi]
                           )
                        << "\n\n        Integral total heat flux [W]      : "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchTotHeatFlux[patchi]
                           )
                        << "\n        Average total heat flux [W/m2]    : "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchTotHeatFlux[patchi]
                           )
                          /gSum(mesh.magSf().boundaryField()[patchi])
                        << nl << endl;
                }
           }
           Info << endl;
        }

    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
