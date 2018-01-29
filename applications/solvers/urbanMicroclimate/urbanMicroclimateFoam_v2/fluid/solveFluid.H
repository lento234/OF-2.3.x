simpleControlFluid simpleControlFluid(fluidRegions[i]);

//simplifiedVegetationModel vegetation(U, thermo.T(), w, Tl); // Vegetation model
//soilVegetationModel vegetation(U, thermo.T(), w, Tl); // Vegetation model

label fluidIteration = 0;

while (simpleControlFluid.loop())
{
    //  Pressure-velocity SIMPLE corrector
    fluidIteration++;
    Info << "Internal iteration number: " << fluidIteration << endl;

    p_rgh.storePrevIter();
    rho.storePrevIter();
    {
        #include "UEqn.H"
        #include "EEqn.H"
        #include "pEqn.H"
        #include "wEqn.H"
    }

    // Solve vegetation energy balance // not really needed every iteration
    if ((fluidIteration % 5) == 0)
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

    // if ((fluidIteration % 50) == 0)
    //     runTime.write();
        
    Info << nl;

}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration << " iterations" << nl << endl;

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();
