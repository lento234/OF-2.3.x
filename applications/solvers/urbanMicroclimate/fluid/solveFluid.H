simpleControlFluid simpleControlFluid(fluidRegions[i]);

label fluidIteration = 1;

while ((simpleControlFluid.loop() || fluidIteration <= minFluidIteration) && fluidIteration <= maxFluidIteration)
{
    //  Pressure-velocity SIMPLE corrector
    Info << "Internal iteration number: " << fluidIteration << endl;

    p_rgh.storePrevIter();
    rho.storePrevIter();
    {
        #include "UEqn.H"
        #include "EEqn.H"
        #include "pEqn.H"
        #include "wEqn.H"
    }

    // Solve vegetation energy balance // every nth iteration
    if ((fluidIteration % 5) == 0) // TODO: solution control -> vegetation solve frequency
    {
        label leafzoneID = -1;
        label rootzoneID = -1;
        forAll(vegetationLeafNames, j)
        {
            if (fluidRegions[i].name() == vegetationLeafNames[j])
                leafzoneID = j;
        }
        forAll(solidRegions, k)
        {
            if (solidRegions[k].name() == vegetationRootNames[leafzoneID])
                rootzoneID = k;
        }
        vegetation.solve(U, thermo, w, wsSolid[rootzoneID]);
    }

    turb.correct();
        
    Info << nl;
    fluidIteration++;

}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration-1 << " iterations" << nl << endl;

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();
