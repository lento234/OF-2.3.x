
const polyBoundaryMesh& patches = mesh.boundaryMesh();
const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();
// const volScalarField::GeometricBoundaryField& Qrb = Qr.boundaryField();
//const volScalarField::GeometricBoundaryField& LADb = LAD.boundaryField();

label nCoarseFaces = 0;      //total number of coarse faces
label nCoarseFacesAll = 0;   //Also includes non-wall faces with greyDiffusive boundary
label nFineFacesTotal = 0;        //total number of fine faces including non-fixedValueFvPatchScalarField patches following advise from bug report

// // Number of boundary faces
// labelList bndFaces(mesh.nFaces()-mesh.nInternalFaces());
//
// label bndI = 0;
// forAll(patches, patchI)
// {
//     const polyPatch& pp = patches[patchI];
//
//     //if (!pp.coupled())
//     if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
//     {
//         forAll(pp, i)
//         {
//             bndFaces[bndI++] = pp.start() + i;
//         }
//     }
// }
// bndFaces.setSize(bndI);

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

// //forAll(LADb, patchI)
// forAll(Qrb, patchI)
// {
//     const polyPatch& pp = patches[patchI];
//     //const fvPatchScalarField& LADpI = LADb[patchI];
//     const fvPatchScalarField& QrpI = Qrb[patchI];
//
//     //if ((isA<fixedValueFvPatchScalarField>(QrpI)) && (pp.size() > 0))
// 	  if ((isA<wallFvPatch>(mesh.boundary()[patchI])) && (pp.size() > 0))
//     {
//         //viewFactorsPatches[count] = QrpI.patch().index();
//         nCoarseFaces += coarsePatches[patchI].size();
// 		    nCoarseFacesAll += coarsePatches[patchI].size();
//         //nFineFaces += patches[patchI].size();
// 		    //count ++;
//         //
// 		    //howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchI].size();
//         //
//         //label i = 0;
//         //for (; i < howManyCoarseFacesPerPatch[countAll]; i++)
//         //{
//         //    sunskyMap_.append(countForMapping);
//         //    countForMapping ++;
//         //}
//          nFineFacesTotal += patches[patchI].size();
//     }
//   	//else if ((isA<fixedValueFvPatchScalarField>(LADpI)) && (pp.size() > 0))
//     else if ((isA<fixedValueFvPatchScalarField>(QrpI)) && (pp.size() > 0))
// 	  {
// 	      nCoarseFacesAll += coarsePatches[patchI].size();
//         //howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchI].size();
//         //
//         //label i = 0;
//         //for (; i < howManyCoarseFacesPerPatch[countAll]; i++)
//         //{
//         //    sunskyMap_.append(countForMapping);
//         //    countForMapping ++;
//         //}
//         nFineFacesTotal += patches[patchI].size();
//    }
//    else
//    {
//        //howManyCoarseFacesPerPatch[countAll] = 0;
//        nFineFacesTotal += patches[patchI].size();
//    }
//    //countAll ++;
// }

Info << "??????????????????????? nFineFacesTotal " << nFineFacesTotal << endl;
Info << "??????????????????????? nCoarseFacesAll " << nCoarseFacesAll << endl;
Info << "??????????????????????? nCoarseFaces " << nCoarseFaces << endl;


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

//surfacesMesh.searchableSurface::write();

surfacesMesh.setField(triSurfaceToAgglom);
