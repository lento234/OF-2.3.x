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
    humidPorousSimpleFoam

Description
    Non-Uniform Porous SimpleFoam with Mass transfer for humidity
        Steady-state solver for incompressible, turbulent flow
        past a non-uniform porous media with transpiration.

Implementation
    Name: Lento Manickathan, <manickathan@arch.ethz.ch>
    Date: Aug, 2015
    Affiliation: EMPA, ETHZ
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H" // added

#include "vegetationModel.H"  // vegetation model added

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    // #include "readGravitationalAcceleration.H" // buoyancy term

    simpleControl simple(mesh);

    vegetationModel vegetation(U, LAD, LAI, T); // Vegetation model

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    // vegetation.testing(U, T, q);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector with decoupled temperature eq.
        {

            #include "UEqn.H" // momentum transport eqn
            #include "TEqn.H" // temperature transport eqn
            #include "pEqn.H" // divergence free velocity field
        }

        // Solve for turbulence
        turbulence->correct();

        // Update leaf temperature
        vegetation.update_Tl(U, T, q);

        // Solve passive scalar transport
        {
            #include "qEqn.H" // specific humidity transport eqn
            // #include "cEqn.H" // CO2 concentration transport eqn
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
