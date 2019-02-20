Random rndGen(653213);

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox
    (
        boundBox(coarseMesh.points(), false)
    ).extend(rndGen, 1e-3)
);

// Dummy bounds dictionary
dictionary dict;
dict.add("bounds", meshBb);
dict.add
(
    "distributionType",
    distributedTriSurfaceMesh::distributionTypeNames_
    [
        distributedTriSurfaceMesh::FROZEN
    ]
);
dict.add("mergeDistance", SMALL);

///////////////////includeWallPatches//////////////////////
labelHashSet includeWallPatches;
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
//    if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
//	if (isA<wallPolyPatch>(pp) || isA<cyclicPolyPatch>(pp)) //ayk
	if (isA<wallPolyPatch>(pp)) //ayk
    {
        includeWallPatches.insert(patchI);
    }
}

labelList triSurfaceToAgglom(5*nFineFacesTotal);

const triSurface localSurface = triangulate
(
    patches,
    includeWallPatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches
);

distributedTriSurfaceMesh surfacesMesh
(
    IOobject
    (
        "wallSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    localSurface,
    dict
);

triSurfaceToAgglom.resize(surfacesMesh.size());

surfacesMesh.searchableSurface::write();

surfacesMesh.setField(triSurfaceToAgglom);
/////////////////////////////////////////////////////////////

///////////////////includeCyclicPatches//////////////////////
labelHashSet includeCyclicPatches;
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
////    if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
    if (isA<wallPolyPatch>(pp) || isA<cyclicPolyPatch>(pp) || isA<processorCyclicPolyPatch>(pp) || isA<emptyPolyPatch>(pp) ) //ayk
//    if (isA<cyclicPolyPatch>(pp)) //ayk    
    {
        includeCyclicPatches.insert(patchI); 
    }
}

labelList triSurfaceToAgglomCyclic(5*nFineFacesTotal);

const triSurface localSurfaceCyclic = triangulate
(
    patches,
    includeCyclicPatches,
    finalAgglom,
    triSurfaceToAgglomCyclic,
    globalNumbering,
    coarsePatches
);

distributedTriSurfaceMesh cyclicMesh
(
    IOobject
    (
        "cyclicSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    localSurfaceCyclic,
    dict
);

triSurfaceToAgglomCyclic.resize(cyclicMesh.size());

cyclicMesh.searchableSurface::write();

cyclicMesh.setField(triSurfaceToAgglomCyclic);
///////////////////////////////////////////////////////////
