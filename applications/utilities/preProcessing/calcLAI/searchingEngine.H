
const polyBoundaryMesh& patches = mesh.boundaryMesh();
const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

label nCoarseFaces = 0;      //total number of coarse faces
label nCoarseFacesAll = 0;   //Also includes non-wall faces with greyDiffusive boundary
label nFineFacesTotal = 0;        //total number of fine faces including non-fixedValueFvPatchScalarField patches following advise from bug report

const volScalarField::GeometricBoundaryField& LADb = LAD.boundaryField();

forAll(LADb, patchI)
{
    const polyPatch& pp = patches[patchI];

    if ((isA<wallFvPatch>(mesh.boundary()[patchI])) && (pp.size() > 0))
    {
      nCoarseFaces += coarsePatches[patchI].size();
      nCoarseFacesAll += coarsePatches[patchI].size();
      nFineFacesTotal += patches[patchI].size();
    }
    else if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
    {
        nCoarseFacesAll += coarsePatches[patchI].size();
        nFineFacesTotal += patches[patchI].size();
    }

}

globalIndex globalNumbering(nCoarseFaces);

Random rndGen(653213);

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox
    (
        boundBox(mesh.points(), false)
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

labelHashSet includePatches;
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
    if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
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

surfacesMesh.setField(triSurfaceToAgglom);
