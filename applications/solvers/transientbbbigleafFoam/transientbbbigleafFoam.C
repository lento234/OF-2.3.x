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
    transientbbbigleafFoam

Description
    "Transient" steady-state solver for buoyant turbulent flow of incompressible
    fluids with vegetation as "big-leaf"

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "fvIOoptionList.H"
#include "transientSimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "transientVegetationModel.H"  // vegetation model added

#include <ctime>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    clock_t tstart;

    transientVegetationModel vegetation(U, T, Yv); // Vegetation model

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        transientSimpleControl transientSimple(mesh);

        int steadyStateIter = 0;

        while (transientSimple.loop())
        {
            tstart = std::clock(); // profiler

            Info << "\nTime = " << runTime.timeName()
                 << ", Steady-state Iteration = " << ++steadyStateIter << endl;

            // Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "TEqn.H"
                #include "pEqn.H"
            }

            // Solve passive scalar transport
            {
                #include "YvEqn.H" // specific humidity transport eqn
            }

            // Solve vegetation energy balance
            vegetation.solve(U, T, Yv);

            // solve k, epsilon
            turbulence->correct();

            // check solution during runtime (will overwrite)
            if ((steadyStateIter % 1000) == 0)
                runTime.write();

            Info << "It took "<< (std::clock()-tstart) / (double)CLOCKS_PER_SEC
                 << " second(s)."<< endl; // profiler
        }

        Info << "\nSteady-state solution converged in "
             << steadyStateIter << " iterations." << endl;

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
