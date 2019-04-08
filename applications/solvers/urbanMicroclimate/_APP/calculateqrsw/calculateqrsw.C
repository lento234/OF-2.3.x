/*---------------------------------------------------------------------------* \
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

Application
    calculateqrsw

Description
    calculateqrsw by Lento Manickathan

Versions
    April 5th

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "fvc.H"
#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMesh.H"
#include "cyclicAMIPolyPatch.H"
#include "triSurfaceTools.H"
#include "mapDistribute.H"
#include "OFstream.H"
#include "meshTools.H"
#include "meshSearch.H"
#include "plane.H"
#include "uindirectPrimitivePatch.H"
#include "DynamicField.H"
#include "IFstream.H"
#include "unitConversion.H"

#include "mathematicalConstants.H"
#include "scalarMatrices.H"
#include "CompactListList.H"
#include "labelIOList.H"
#include "labelListIOList.H"
#include "scalarListIOList.H"
#include "vectorListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"

#include "singleCellFvMesh.H"
#include "interpolation.H"
#include "IOdictionary.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
#include "unitConversion.H"

using namespace Foam;

// calculate the end point for a ray hit check
point calcEndPoint
(
    const point &start,
    const point &n2,
    const point &pminO,
    const point &pmaxO
)
{
  scalar ix = 0; scalar iy = 0; scalar iz = 0;

  if (n2.x() > 0.0)
    ix = (pmaxO.x() - start.x())/n2.x();
  else if (n2.x() < 0.0)
    ix = (pminO.x() - start.x())/n2.x();
  else
    ix = VGREAT;

  if (n2.y() > 0.0)
    iy = (pmaxO.y() - start.y())/n2.y();
  else if (n2.y() < 0.0)
    iy = (pminO.y() - start.y())/n2.y();
  else iy = VGREAT;

  if (n2.z() > 0.0)
    iz = (pmaxO.z() - start.z())/n2.z();
  else if (n2.z() < 0.0)
    iz = (pminO.z() - start.z())/n2.z();
  else
    iz = VGREAT;

  // closest edg direction
  scalar i = min(ix, min(iy, iz));

  return 0.9999*i*n2 + start;
}

triSurface triangulate
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& includePatches,
    const labelListIOList& finalAgglom,
    labelList& triSurfaceToAgglom,
    const globalIndex& globalNumbering,
    const polyBoundaryMesh& coarsePatches
)
{
    const polyMesh& mesh = bMesh.mesh();

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles
    (
        mesh.nFaces() - mesh.nInternalFaces()
    );

    label newPatchI = 0;
    label localTriFaceI = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchI = iter.key();
        const polyPatch& patch = bMesh[patchI];
        const pointField& points = patch.points();

        label nTriTotal = 0;

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            faceList triFaces(f.nTriangles(points));

            label nTri = 0;

            f.triangles(points, nTri, triFaces);

            forAll(triFaces, triFaceI)
            {
                const face& f = triFaces[triFaceI];

                triangles.append(labelledTri(f[0], f[1], f[2], newPatchI));

                nTriTotal++;

                triSurfaceToAgglom[localTriFaceI++] = globalNumbering.toGlobal
                (
                    Pstream::myProcNo(),
                    finalAgglom[patchI][patchFaceI]
                  + coarsePatches[patchI].start()
                );
            }
        }

        newPatchI++;
    }

    triSurfaceToAgglom.resize(localTriFaceI-1);

    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().setSize(newPatchI);

    newPatchI = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchI = iter.key();
        const polyPatch& patch = bMesh[patchI];

        surface.patches()[newPatchI].index() = patchI;
        surface.patches()[newPatchI].name() = patch.name();
        surface.patches()[newPatchI].geometricType() = patch.type();

        newPatchI++;
    }

    return surface;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // start timer
    clock_t tstart = std::clock();
    clock_t tstartStep;

    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    Info << "timeDirs: " << timeDirs << endl;
    runTime.setTime(timeDirs[0], 0);

    #include "createNamedMesh.H"

    volScalarField Qr
    (
      IOobject
      (
          "Qr",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      ),
      mesh
    );

    //wordList boundaryTypes = Qr.boundaryField().types();
    wordList boundaryTypes(Qr.boundaryField().types().size(), "zeroGradient");

    Info << "boundaryTypes = " << boundaryTypes << endl;

    // Read sunPosVector list
    vectorIOList sunPosVector
    (
       IOobject
       (
            "sunPosVector",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );

    scalarIOList IDN // direct solar radiation intensity flux
    (
       IOobject
       (
            "IDN",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );

    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    singleCellFvMesh coarseMesh
    (
        IOobject
        (
            mesh.name(),
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        finalAgglom
    );

    #include "readGravitationalAcceleration.H"
    Info << "Gravity is = " << g << endl;

    const vector ez = - g.value()/mag(g.value());
    Info << "Vertical vector : " << ez << endl;

    // Mesh setup
    int nMeshCells = mesh.cells().size();

    // Mesh cell centers
    pointField pmeshC = mesh.C();

    // mesh bounding box
    point pminO = gMin(pmeshC);
    point pmaxO = gMax(pmeshC);

    // Set up searching engine for obstacles (copied from Aytac)
    #include "searchingEngine.H"

    DynamicList<label> rayStartFace(nCoarseFaces + 0.01*nCoarseFaces);

    tstartStep = std::clock();

    forAll(sunPosVector, vectorID)
    {

        Info << nl << "Time = " << runTime.timeName() << endl;

        // start clock
        tstartStep = std::clock();

        // Volume vector field qrsw 
        volVectorField qrswi
        (
            IOobject
            (
               "qrsw",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE 
            ),
            mesh,
            dimensionedVector("0", dimensionSet(1,0,-2,0,0,0,0), vector::zero),
            boundaryTypes
        );
        
        // sunPosVector i
        vector n2 = sunPosVector[vectorID];
        n2 /= mag(n2);

        // only if sun is above the horizon
        if ( (n2 & ez) > 0)
        {

            DynamicField<point> startListCells(nMeshCells);
            DynamicField<point> endListCells(nMeshCells);
            List<pointIndexHit> pHitListCells(nMeshCells);

            // Calculate solar short-wave radiation vector field
            forAll(qrswi, cellI)
            {
                point starti = pmeshC[cellI];
                point endi = calcEndPoint(starti, n2, pminO, pmaxO);
                startListCells.append(starti);
                endListCells.append(endi);
            }

            surfacesMesh.findLine(startListCells, endListCells, pHitListCells);

            forAll(qrswi, cellI)
            {
                if (!pHitListCells[cellI].hit())
                {
                    qrswi[cellI] = -n2*IDN[vectorID];
                }
            }
        }

        // Export --- temp.
        //Info << "Info: Exporting step " << vectorID << endl;
        
        qrswi.correctBoundaryConditions();
        qrswi.write();

        runTime++;

        
        Info << "Solar ray direction " << vectorID
             << ", It took " ;
        printf("%.3f", (std::clock()-tstartStep) / static_cast<double>(CLOCKS_PER_SEC));
        //printf("%.3f", (std::clock()-tstartStep) / (double)CLOCKS_PER_SEC);
        Info << " second(s)."<< endl;
        

    }

    // Status Info
    Info << "\nTotal time took: "
         //<< (std::clock()-tstart) / (double)CLOCKS_PER_SEC
         << (std::clock()-tstart) / static_cast<double>(CLOCKS_PER_SEC)
         <<" second(s).\n"<< endl;

    Info<< "End\n" << endl;
    return 0;

}


// ************************************************************************* //
