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

///////////////////includePatches//////////////////////
labelHashSet includePatches;
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
    if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp) && !isA<emptyPolyPatch>(pp))
    {
        includePatches.insert(patchI);
    }
}

labelList triSurfaceToAgglom(5*nFineFacesTotal);

const triSurface localSurface = triangulate
(
    patches,
    includePatches,
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
/*labelHashSet includeCyclicPatches;
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
    if (!pp.coupled())   
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

cyclicMesh.setField(triSurfaceToAgglomCyclic);*/
///////////////////////////////////////////////////////////