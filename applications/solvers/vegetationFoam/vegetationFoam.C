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
    vegetationFoam

Description
    Solves for air flow through vegetation
          Written by Lento Manickathan, September 2017, ETH Zurich/Empa

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulenceModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "solarLoadModel.H"             // solar radiation model by aytac
#include "simpleControlFluid.H"         // simple control by aytac
#include "fvIOoptionList.H"
#include "coordinateSystem.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "simplifiedVegetationModel.H"  // vegetation model by Lento added

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"

    #include "createFluidFields.H"
    //#include "createSolidFields.H"

    #include "initContinuityErrs.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        forAll(fluidRegions, i)
        {
            if (fluidRegions[i].name() == "vegetation")
			      {
				          Info<< "\nVegetation region found..." << endl;
				          #include "setVegetationFields.H"

				          Info << "Updating T boundary fields..." << endl;
				          thermo.T().correctBoundaryConditions();

                  Info << "Updating long-wave radiation heat transfer for region: "
                       << fluidRegions[i].name() << endl;
				          rad.correct();

                  Info << "Updating short-wave radiation heat transfer for region: "
                       << fluidRegions[i].name() << endl;
				          sol.correct();
			      }
			      else
			      {
	               Info << "\nSolving for fluid region "
	                    << fluidRegions[i].name() << endl;
	               #include "setRegionFluidFields.H"
	               #include "readFluidMultiRegionSIMPLEControls.H"
	               #include "solveFluid.H"

	          }
        }


        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        //scalar storeFluidDeltaT = runTime.deltaT().value();


        forAll(solidRegions, i)
		    {
			      Info << "\nSolving for solid region "
				         << solidRegions[i].name()
                 << " not implemented !!!!!" << endl;
			      //#include "setRegionSolidFields.H"
			      //#include "readSolidMultiRegionSIMPLEControls.H"
			      //#include "solveSolid.H"
		    }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
