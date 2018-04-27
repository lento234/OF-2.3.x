/*---------------------------------------------------------------------------*\
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

Global
    readFluidControls

Description
    Read the control parameters used in the fluid

\*---------------------------------------------------------------------------*/

label maxFluidIteration =
    runTime.controlDict().lookupOrDefault<label>("maxFluidIteration", 10000);

label minFluidIteration =
    runTime.controlDict().lookupOrDefault<label>("minFluidIteration", 1);

label vegetationSolveFreq =
    runTime.controlDict().lookupOrDefault<label>("vegetationSolveFreq", 10);    

bool debugExport =
    runTime.controlDict().lookupOrDefault<bool>("debugExport", false);

Info << "maxFluidIteration: " << maxFluidIteration << endl;

Info << "vegetationSolveFreq: " << vegetationSolveFreq << endl;

Info << "debugExport: " << debugExport << endl;


// ************************************************************************* //
