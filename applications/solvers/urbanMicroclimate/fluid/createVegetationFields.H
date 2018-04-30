

    // Info
    const wordList vegetationLeafNames(rp["vegetation_leafzone"]);
    const wordList vegetationRootNames(rp["vegetation_rootzone"]);

    if (vegetationLeafNames.size() != vegetationRootNames.size())
    {
        Info << "WARNING:: vegetationLeafNames.size(): " << vegetationLeafNames.size()
             << ", vegetationRootNames.size()" << vegetationRootNames.size()
             << endl;
        return 1;
    }

    label fluidRegionID = -1;
    label solidRegionID = -1;

    label nVegDomains = vegetationLeafNames.size();

    PtrList<soilVegetationModel> vegetationFluid(nVegDomains); // Vegetation model

    forAll(vegetationLeafNames, i)
    {
        forAll(fluidRegions, j)
        {
            if (fluidRegions[j].name() == vegetationLeafNames[i])
                fluidRegionID = j;
        }

        forAll(solidRegions, k)
        {
            if (solidRegions[k].name() == vegetationRootNames[i])
                solidRegionID = k;
        }

        
        Info <<"    Adding vegetation\n" << endl;

        Info << "       fluidRegionID: " << fluidRegionID
             << ", solidRegionID: " << solidRegionID << endl;


        vegetationFluid.set
        (
            i,
            new soilVegetationModel
            (
                UFluid[fluidRegionID],
                thermoFluid[fluidRegionID],
                wFluid[fluidRegionID],
                cFluid[fluidRegionID],
                //TsSolid[solidRegionID],
                //wsSolid[solidRegionID],
                pcSolid[solidRegionID]
            )
        );        

    }