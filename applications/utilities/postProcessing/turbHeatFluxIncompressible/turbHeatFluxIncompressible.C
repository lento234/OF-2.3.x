/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2010-08-02 Eelco van Vliet: 1st public version of wallHeatFluxIncompressible:
  http://www.cfd-online.com/Forums/openfoam-solving/66705-wallheatflux-bc-not-constant-after-restart.html#post269812

2012-05-21 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post362191
  «modified the standard wallHeatflux utility which comes default with OF into
  a version for incompressible flows. Also removed a bug out of the code.»

2012-06-26 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post368330
  «p is now not required anymore.»

2014-06-22: Bruno Santos: Adapted to OpenFOAM 2.2.x.

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
    turbHeatFluxIncompressible

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.
    Based on wallHeatFlux with changes to allow it on incompressible flows
    Also removed a bug at the typeid checkline
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"
        #include "readTransportProperties.H"

         // update the turbulence fields
        turbulence->read();

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
            alphaEff=turbulence->nu()/Pr+alphat;
        }
        else
        {
            Info<< "\nRead effective heat conductivity alphaEff" << endl;
        }

      // Volume fields

      // grad T
      volVectorField gradT("gradT", fvc::grad(T));
      gradT.correctBoundaryConditions();

      // turbulent heat flux  Qturb = alphaeff*rho*Cp*gradT
      volVectorField turbHeatFlux("turbHeatFlux", alphaEff*Cp*rho*gradT);
      //turbHeatFlux.correctBoundaryConditions();

      // surface normal heat flux
      surfaceScalarField turbHeatFluxNormal = fvc::interpolate(turbHeatFlux) & mesh.Sf()/mesh.magSf();

      const surfaceScalarField::GeometricBoundaryField& patchturbHeatFlux =
               turbHeatFluxNormal.boundaryField();

      Info<< "\nWall heat fluxes " << endl;
      forAll(patchturbHeatFlux, patchi)
      {
          if ( (!isA<emptyFvPatch>(mesh.boundary()[patchi])) &&
             (mesh.magSf().boundaryField()[patchi].size() > 0) )
          {
              Info<< mesh.boundary()[patchi].name()
                  << ": Total "
                  << gSum
                     (
                         mesh.magSf().boundaryField()[patchi]
                        *patchturbHeatFlux[patchi]
                     )
                  << " [W] over "
                  << gSum(mesh.magSf().boundaryField()[patchi])
                  << " [m2] ("
                  << gSum
                     (
                         mesh.magSf().boundaryField()[patchi]
                        *patchturbHeatFlux[patchi]
                     )
                    /gSum(mesh.magSf().boundaryField()[patchi])
                  << " [W/m2])"
                  << endl;
        }
    }
    Info << endl;

    // export
    gradT.write();
    turbHeatFlux.write();

    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
