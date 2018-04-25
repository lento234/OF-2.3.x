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

///////////////////stl original//////////////////////
point moveSTL = vector(0,0,0);

const triSurface localSurface = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches,
    moveSTL
);
/////////////////////////////////////////////////////

///////////////////stl XPos//////////////////////
moveSTL = vector(maxList[0][0],0,0);

const triSurface localSurfaceXPos = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches,
    moveSTL
);
/////////////////////////////////////////////////////
                
///////////////////stl XNeg//////////////////////
moveSTL = vector(-maxList[0][0],0,0);

const triSurface localSurfaceXNeg = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches,
    moveSTL
);
/////////////////////////////////////////////////////

//////////////COMBINING step 1////////////////////
    //points
    const pointField& points1 = localSurface.points();
    const pointField& points2 = localSurfaceXPos.points();
    const pointField& points3 = localSurfaceXNeg.points();    
    vectorField pointsAll(points1.size() + points2.size() + points3.size());        
    label pointi = 0;
    // Copy points1 into pointsAll
    forAll(points1, point1i)
    {
        pointsAll[pointi++] = points1[point1i];
    }
    // Add surface2 points
    forAll(points2, point2i)
    {
        pointsAll[pointi++] = points2[point2i];
    }
    // Add surface3 points
    forAll(points3, point3i)
    {
        pointsAll[pointi++] = points3[point3i];
    }
    
    //faces
    // Copy triangles1 into trianglesAll
    List<labelledTri> facesAll(localSurface.size() + localSurfaceXPos.size() + localSurfaceXNeg.size());
    label trianglei = 0;
    forAll(localSurface, facei)
    {
        facesAll[trianglei++] = localSurface[facei];
    }
    label nRegions1 = localSurface.patches().size();
    // Add (renumbered) surface2 triangles
    forAll(localSurfaceXPos, facei)
    {
        const labelledTri& tri = localSurfaceXPos[facei];

        labelledTri& destTri = facesAll[trianglei++];
        destTri[0] = tri[0] + points1.size();
        destTri[1] = tri[1] + points1.size();
        destTri[2] = tri[2] + points1.size();
        destTri.region() = tri.region() + nRegions1;
    }    
    label nRegions2 = localSurfaceXPos.patches().size();
    // Add (renumbered) surface3 triangles
    forAll(localSurfaceXNeg, facei)
    {
        const labelledTri& tri = localSurfaceXNeg[facei];

        labelledTri& destTri = facesAll[trianglei++];
        destTri[0] = tri[0] + points1.size() + points2.size();
        destTri[1] = tri[1] + points1.size() + points2.size();
        destTri[2] = tri[2] + points1.size() + points2.size();
        destTri.region() = tri.region() + nRegions1 + nRegions2;
    }    
    label nRegions3 = localSurfaceXNeg.patches().size();
    
    //patches
    geometricSurfacePatchList newPatches;
    newPatches.setSize(nRegions1 + nRegions2 + nRegions3);
    label newPatchi = 0;    
    forAll(localSurface.patches(), patchi)
    {
        newPatches[newPatchi++] = localSurface.patches()[patchi];
    }        
    forAll(localSurfaceXPos.patches(), patchi)
    {
        newPatches[newPatchi++] = localSurfaceXPos.patches()[patchi];
    }    
    forAll(localSurfaceXNeg.patches(), patchi)
    {
        newPatches[newPatchi++] = localSurfaceXNeg.patches()[patchi];
    }    
    
    triSurface combinedSurf;
    combinedSurf = triSurface(facesAll, newPatches, pointsAll);

    // Merge all common points and do some checks
    //combinedSurf.cleanup(true);

    Info<< "Merged surface:" << endl;

    combinedSurf.writeStats(Info);
/////////////////////////////////////////////////////        

//////////////TRANSLATING////////////////////    
    pointField points1_(combinedSurf.points());
    pointField points2_(combinedSurf.points());
    vector v1 = vector(0,0,maxList[0][2]);
    vector v2 = vector(0,0,-maxList[0][2]);
    points1_ += v1; 
    points2_ += v2; 
    
    triSurface combinedSurf2 = combinedSurf;
    combinedSurf2.movePoints(points1_);
    Info<< "Translated surface:" << endl;

    combinedSurf2.writeStats(Info);
    
    triSurface combinedSurf3 = combinedSurf;
    combinedSurf3.movePoints(points2_);
    Info<< "Translated surface:" << endl;

    combinedSurf2.writeStats(Info);    
/////////////////////////////////////////////////////    

//////////////COMBINING step 2////////////////////
    //points
    const pointField& points1__ = combinedSurf.points();
    const pointField& points2__ = combinedSurf2.points();
    const pointField& points3__ = combinedSurf3.points();    
    vectorField pointsAll__(points1__.size() + points2__.size() + points3__.size());        
    pointi = 0;
    // Copy points1 into pointsAll
    forAll(points1__, point1i)
    {
        pointsAll__[pointi++] = points1__[point1i];
    }
    // Add surface2 points
    forAll(points2__, point2i)
    {
        pointsAll__[pointi++] = points2__[point2i];
    }
    // Add surface3 points
    forAll(points3__, point3i)
    {
        pointsAll__[pointi++] = points3__[point3i];
    }
    
    //faces
    // Copy triangles1 into trianglesAll
    List<labelledTri> facesAll__(combinedSurf.size() + combinedSurf2.size() + combinedSurf3.size());
    trianglei = 0;
    forAll(combinedSurf, facei)
    {
        facesAll__[trianglei++] = combinedSurf[facei];
    }
    nRegions1 = combinedSurf.patches().size();
    // Add (renumbered) surface2 triangles
    forAll(combinedSurf2, facei)
    {
        const labelledTri& tri = combinedSurf2[facei];

        labelledTri& destTri = facesAll__[trianglei++];
        destTri[0] = tri[0] + points1__.size();
        destTri[1] = tri[1] + points1__.size();
        destTri[2] = tri[2] + points1__.size();
        destTri.region() = tri.region() + nRegions1;
    }    
    nRegions2 = combinedSurf2.patches().size();
    // Add (renumbered) surface3 triangles
    forAll(combinedSurf3, facei)
    {
        const labelledTri& tri = combinedSurf3[facei];

        labelledTri& destTri = facesAll__[trianglei++];
        destTri[0] = tri[0] + points1__.size() + points2__.size();
        destTri[1] = tri[1] + points1__.size() + points2__.size();
        destTri[2] = tri[2] + points1__.size() + points2__.size();
        destTri.region() = tri.region() + nRegions1 + nRegions2;
    }    
    nRegions3 = combinedSurf3.patches().size();
    
    //patches
    geometricSurfacePatchList newPatches__;
    newPatches__.setSize(nRegions1 + nRegions2 + nRegions3);
    newPatchi = 0;    
    forAll(combinedSurf.patches(), patchi)
    {
        newPatches__[newPatchi++] = combinedSurf.patches()[patchi];
    }        
    forAll(combinedSurf2.patches(), patchi)
    {
        newPatches__[newPatchi++] = combinedSurf2.patches()[patchi];
    }    
    forAll(combinedSurf3.patches(), patchi)
    {
        newPatches__[newPatchi++] = combinedSurf3.patches()[patchi];
    }    
    
    triSurface combinedSurfFinal;
    combinedSurfFinal = triSurface(facesAll__, newPatches__, pointsAll__);

    // Merge all common points and do some checks
    //combinedSurf.cleanup(true);

    Info<< "Merged Final surface:" << endl;

    combinedSurfFinal.writeStats(Info);
/////////////////////////////////////////////////////    

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
    combinedSurfFinal,
    dict
);

triSurfaceToAgglom.resize(surfacesMesh.size());

surfacesMesh.searchableSurface::write();

surfacesMesh.setField(triSurfaceToAgglom);