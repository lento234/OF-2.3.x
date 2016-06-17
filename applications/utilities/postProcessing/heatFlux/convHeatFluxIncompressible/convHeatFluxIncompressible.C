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
    convHeatFluxIncompressible

Description
    Calculates and writes the convective heat flux for fluid domain and all
    patches as the boundary field.
    Prints the integrated flux for all (non-empty) patches.
    Based on wallHeatFlux with changes to allow it on incompressible flows
    Also removed a bug at the typeid checkline
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    argList::addBoolOption
    (
      "print",
      "print turbulent heat fluxes at the boundaries"
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

        // calculate convective heat flux Qconv
        volVectorField convHeatFlux("convHeatFlux", U*(T-TRef)*Cp*rho);

        Info<< "\nWriting convHeatFlux field" << nl << endl;
        convHeatFlux.write();

        // print results
        if (printResults)
        {
            surfaceScalarField convHeatFluxNormal = fvc::interpolate(convHeatFlux) & mesh.Sf()/mesh.magSf();

            const surfaceScalarField::GeometricBoundaryField& patchConvHeatFlux =
                convHeatFluxNormal.boundaryField();

            Info<< "\nConvective heat fluxes at the boundaries" << endl;
            forAll(patchConvHeatFlux, patchi)
            {
                if ( (!isA<emptyFvPatch>(mesh.boundary()[patchi])) &&
                     (mesh.boundary()[patchi].size() > 0) )
                {
                    Info<< "    "
                        << mesh.boundary()[patchi].name()
                        << "\n        Total heat flux [W]  : "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchConvHeatFlux[patchi]
                           )
                        << "\n        Total area [m2]      : "
                        << gSum(mesh.magSf().boundaryField()[patchi])
                        << "\n        Avg. heat flux [W/m2]: "
                        << gSum
                           (
                               mesh.magSf().boundaryField()[patchi]
                              *patchConvHeatFlux[patchi]
                           )
                          /gSum(mesh.magSf().boundaryField()[patchi])
                        << endl;
                }
           }
           Info << endl;
        }

    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
