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
    calcLAI

Description
    calcLAI by Lento Manickathan

Versions
    May   - v1
    July  - v2
    Aug   - v3, v4

\*---------------------------------------------------------------------------*/

#include <sstream>
#include <string>
#include <iostream>

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
//#include "regionProperties.H"
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
//#include "treeDataFace.H"
#include "unitConversion.H"
//#include "fvIOoptionList.H"

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


// trilinear interpolation
scalar interp3D
(
    const point &ptemp,
    const pointField &pInterp,
    const scalarField &valInterp,
    const int &nx,
    const int &ny
) //point &dp,
{

  // point &dp, int &i0, int &j0, int &k0, int &nx, int &ny, int &nz)
  point pmin = pInterp[0];
  point dp   = pInterp[(nx*ny) + nx + 1] - pInterp[0];

  // Offset from p min.
  double xp = ptemp.x()-pmin.x();
  double yp = ptemp.y()-pmin.y();
  double zp = ptemp.z()-pmin.z();

  // Determine index of lower bound
  int i0 = floor(xp/dp.x());
  int j0 = floor(yp/dp.y());
  int k0 = max(0, floor(zp/dp.z())); // if (k0 < 0) k0 = 0;


  // indices
  int i000 = (nx*ny)*k0 + j0*nx + i0;
  int i100 = (nx*ny)*k0 + j0*nx + i0+1;
  int i010 = (nx*ny)*k0 + (j0+1)*nx + i0;
  int i110 = (nx*ny)*k0 + (j0+1)*nx + i0+1;
  int i001 = (nx*ny)*(k0+1) + j0*nx + i0;
  int i101 = (nx*ny)*(k0+1) + j0*nx + i0+1;
  int i011 = (nx*ny)*(k0+1) + (j0+1)*nx + i0;
  int i111 = (nx*ny)*(k0+1) + (j0+1)*nx + i0+1;

  point pInterp000 = pInterp[i000];

  // Bilinear interpolation
  double xd = (ptemp.x()-pInterp000.x())/dp.x();
  double yd = (ptemp.y()-pInterp000.y())/dp.y();
  double zd = (ptemp.z()-pInterp000.z())/dp.z();

  // Interpolation in x-dir
  double c00 = valInterp[i000]*(1-xd) + valInterp[i100]*xd;
  double c01 = valInterp[i001]*(1-xd) + valInterp[i101]*xd;
  double c10 = valInterp[i010]*(1-xd) + valInterp[i110]*xd;
  double c11 = valInterp[i011]*(1-xd) + valInterp[i111]*xd;

  // Interpolation in y-dir
  double c0 = c00*(1.0-yd) + c10*yd;
  double c1 = c01*(1.0-yd) + c11*yd;

  // Interpolate in z-dir
  return c0*(1.0-zd) + c1*zd;

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

void calcVegBBOX
(
    const pointField &pmeshC,
    const volScalarField &LAD,
    point &pmin,
    point &pmax,
    int &vegetationCell
)
{
    point ptemp;

    forAll(LAD, cellI)
    {
        // where vegetation is present
        if (LAD[cellI] > 10*SMALL)
        {
            ptemp = pmeshC[cellI];
            pmin = min(pmin,ptemp);
            pmax = max(pmax,ptemp);
            vegetationCell = cellI;
        }
    }

    List<point> pmin_(Pstream::nProcs());
    List<point> pmax_(Pstream::nProcs());

    pmin_[Pstream::myProcNo()] = pmin;
    pmax_[Pstream::myProcNo()] = pmax;
    Pstream::gatherList(pmin_);
    Pstream::scatterList(pmin_);
    Pstream::gatherList(pmax_);
    Pstream::scatterList(pmax_);

    pmin = gMin(pmin_);
    pmax = gMax(pmax_);

}

void interpfvMeshToCartesian
(
    const fvMesh& mesh,
    const volScalarField &LAD,
    point &pmin,
    point &pmax,
    const int &vegetationCell,
    pointField &pInterp,
    scalarField &LADInterp,
    int &nx,
    int &ny,
    int &nz,
    point &dp
)
{

    // Define search mesh
    meshSearch ms(mesh);

    // Define LAD interpolator to arbitrary locations
    dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes"); // Read interpolation scheme
    autoPtr<interpolation<scalar> > LAD_interpolator = interpolation<scalar>::New(interpolationDict, LAD); // Define interpolator

    // Define cartesian interpolation grid
    scalar minCellV = gMin(mesh.V()); // Cartesian mesh resolution (determine from minimum cell size)
    scalar minCellL = Foam::pow(minCellV, 1.0/3.0);

    dp = vector(minCellL,minCellL,minCellL); // grid spacing

    // Extend the cartesian grid to include vegetation
    pmin -= 5*dp;
    pmax += 5*dp;

    // Define cartesian grid size
    nx = ceil( (pmax.x()-pmin.x()) / dp.x()) + 1;
    ny = ceil( (pmax.y()-pmin.y()) / dp.y()) + 1;
    nz = ceil( (pmax.z()-pmin.z()) / dp.z()) + 1;

    // Generate cartesian interpolation grid
    // coordinates
    pInterp.setSize(nx*ny*nz, point::zero);
    LADInterp.setSize(nx*ny*nz, pTraits<scalar>::zero);

    ////////////////////////////////////////////////////////////////////////////

    /////////////// (Step 2) Interpolate LAD onto cartesian interpolation mesh

    int cellIndex;
    int pIndex;
    point ptemp;

    for (int k=0; k < nz; k++)
    {
      for (int j=0; j < ny; j++)
      {
        for (int i=0; i < nx; i++)
        {
          // index of node p
          pIndex = (nx*ny)*k + j*nx + i;

          // x,y,z coordinates in rotated coordinate system
          pInterp[pIndex].x() = pmin.x() + i*dp.x();
          pInterp[pIndex].y() = pmin.y() + j*dp.y();
          pInterp[pIndex].z() = pmin.z() + k*dp.z();

          // coordinate of point in original coordinate system
          ptemp = pInterp[pIndex]; //ptemp = transform(Tinv, pInterp[pIndex]);

          // Find intersecting cell
          cellIndex = ms.findCell(ptemp,vegetationCell,true); // faster

          // if point is inside domain
          if (cellIndex != -1)
              LADInterp[pIndex] = LAD[cellIndex]; //LAD_interpolator->interpolate(pInterp[pIndex],cellIndex); //TODO: some bug with interpolator

        }
      }
    }

    List<scalarField> LADInterp_(Pstream::nProcs());
    LADInterp_[Pstream::myProcNo()] = LADInterp;
    Pstream::gatherList(LADInterp_);

    if (Pstream::master())
    {
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI > 0)
            {
                LADInterp_[Pstream::myProcNo()] += LADInterp_[procI];
            }
        }
        LADInterp = LADInterp_[Pstream::myProcNo()];
    }

    Pstream::scatter(LADInterp);

    // update maximum point
    pmax = gMax(pInterp);

}

void interpcartesianToRotCartesian
(
    point &pminRot,
    point &pmaxRot,
    pointField &pInterpRot,
    scalarField &LADInterpRot,
    int &nxRot,
    int &nyRot,
    int &nzRot,
    const tensor &Tinv,
    const pointField &pInterp,
    const scalarField &LADInterp,
    const point& pmin,
    const point& pmax,
    const int &nx,
    const int &ny,
    const point &dp
)
{
    pminRot -= 5*dp;
    pmaxRot += 5*dp;

    // Define rotated cartesian grid size
    nxRot = ceil( (pmaxRot.x()-pminRot.x()) / dp.x()) + 1;
    nyRot = ceil( (pmaxRot.y()-pminRot.y()) / dp.y()) + 1;
    nzRot = ceil( (pmaxRot.z()-pminRot.z()) / dp.z()) + 1;

    // Generate rotated cartesian interpolation grid
    // coordinates
    pInterpRot.setSize(nxRot*nyRot*nzRot, point::zero);
    LADInterpRot.setSize(nxRot*nyRot*nzRot, pTraits<scalar>::zero);

    ////////////////////////////////////////////////////////////////////
    // Interpolate onto rotated cartesian grid
    int pIndex;
    point ptemp;

    for (int k=0; k < nzRot; k++)
    {
      for (int j=0; j < nyRot; j++)
      {
        for (int i=0; i < nxRot; i++)
        {
          // index of node p
          pIndex = (nxRot*nyRot)*k + j*nxRot + i;

          // x,y,z coordinates in rotated coordinate system
          pInterpRot[pIndex].x() = pminRot.x() + i*dp.x();
          pInterpRot[pIndex].y() = pminRot.y() + j*dp.y();
          pInterpRot[pIndex].z() = pminRot.z() + k*dp.z();

          // coordinate of point in original coordinate system
          ptemp = transform(Tinv, pInterpRot[pIndex]);

          // If point is within the bbox of original cartesian grid
          if ( (ptemp.x() >= pmin.x()) && (ptemp.x() <= pmax.x()) &&
               (ptemp.y() >= pmin.y()) && (ptemp.y() <= pmax.y()) &&
               (ptemp.z() >= pmin.z()) && (ptemp.z() <= pmax.z()) )
          LADInterpRot[pIndex] = interp3D(ptemp, pInterp, LADInterp, nx, ny);

        }
      }
    }

    // update maximum point
    pmaxRot = gMax(pInterpRot);

}

void integrateLAD
(
    const scalarField &LADInterpRot,
    const int &nxRot,
    const int &nyRot,
    const int &nzRot,
    const point &dp,
    scalarField &LAIInterpRot
)
{
    int pIndex, pIndexkp1;
    LAIInterpRot.setSize(nxRot*nyRot*nzRot, pTraits<scalar>::zero);

    for (int i=0; i < nxRot; i++)
    {
      for (int j=0; j < nyRot; j++)
      {
        for (int k=(nzRot-2); k>=0; k--)
        {
          // lower and upper row index
          pIndex = (nxRot*nyRot)*k + j*nxRot + i;
          pIndexkp1 = (nxRot*nyRot)*(k+1) + j*nxRot + i;
          // trapezoidal integration
          LAIInterpRot[pIndex] = LAIInterpRot[pIndexkp1] + 0.5*(LADInterpRot[pIndex]+LADInterpRot[pIndexkp1])*dp.z();

        }
      }
    }

}

/*
// calculate divergence using finite difference
void calcDiv
(
  const scalarField &qrswInterpRot,
  scalarField &divqrswInterpRot,
  const int &nxRot,
  const int &nyRot,
  const int &nzRot,
  const point &dp
)
{
    int pIndex, pIndexkp1, pIndexkm1;
    divqrswInterpRot.setSize(nxRot*nyRot*nzRot, pTraits<scalar>::zero);

    for (int k=0; k < nzRot; k++)
    {
      for (int j=0; j < nyRot; j++)
      {
        for (int i=0; i < nxRot; i++)
        {
            // p indicies
            pIndex = (nxRot*nyRot)*k + j*nxRot + i; // k
            pIndexkp1 = (nxRot*nyRot)*(k+1) + j*nxRot + i; // p+1 index (forward)
            pIndexkm1 = (nxRot*nyRot)*(k-1) + j*nxRot + i; // p-1 index (backward)

            // forward (k==0), backward (k==end), central (0<k<end)
            if (k == 0)
                divqrswInterpRot[pIndex] = (qrswInterpRot[pIndexkp1] - qrswInterpRot[pIndex])/dp.z();
            else if (k==(nzRot-1))
                divqrswInterpRot[pIndex] = (qrswInterpRot[pIndex] - qrswInterpRot[pIndexkm1])/dp.z();
            else
                divqrswInterpRot[pIndex] = (qrswInterpRot[pIndexkp1] - 2.0*qrswInterpRot[pIndex] + qrswInterpRot[pIndexkm1])/(dp.z()*dp.z());

            //Info << "k = " << k << ", divqrsw = " << divqrswInterpRot[pIndex] << endl;

        }
      }
    }

    Info << "divqrswInterpRot: " << divqrswInterpRot.size() << endl;
    Info << "gMax(divqrswInterpRot): " << gMax(divqrswInterpRot) << endl;
    Info << "gMin(divqrswInterpRot): " << gMin(divqrswInterpRot) << endl;

}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    // start timer
    clock_t tstart = std::clock();
    //clock_t tstartlocal = std::clock();
    clock_t tstartStep;

    timeSelector::addOptions();
    #include "addRegionOption.H"

    Foam::argList::addOption
    (
         "writeFields",
         "",
         "write volScalarFields of all time steps"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // instantList timeDirs = timeSelector::select0(runTime, args);
    //
    // Info << "timeDirs: " << timeDirs << endl;
    //runTime.setTime(label(172800), 0);

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

    wordList boundaryTypes = Qr.boundaryField().types();

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

    vectorListIOList qrswList
    (
        IOobject
        (
            "qrsw",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector.size()
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

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Global setup
    Info << "Creating interpolation setup...\n" << endl;

    // Mesh setup
    int nMeshCells = mesh.cells().size();

    vector n1(0,0,1); // original vector
    n1 /= mag(n1);

    // Mesh cell centers
    pointField pmeshC = mesh.C();

    // mesh bounding box
    point pminO = gMin(pmeshC);
    point pmaxO = gMax(pmeshC);


    // Set up searching engine for obstacles (copied from Aytac)
    #include "searchingEngine.H"

    DynamicList<label> rayStartFace(nCoarseFaces + 0.01*nCoarseFaces);

    // Generate dummy data
    scalarList zeroList_nMeshCells(nMeshCells, 0.0);
    vectorList zeroListVectors_nMeshCells(nMeshCells, point::zero);


    ////// Define LAD interpolator to arbitrary locations
    tstartStep = std::clock();

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    scalarList zeroList_nCoarseFacesAll(nCoarseFacesAll, 0.0);


    ////////////////////////////////////////////////////////////////////////////
    // Status Info
    // Info << "Info: Time, Initialization took: "
    //      << (std::clock()-tstart) / (double)CLOCKS_PER_SEC
    //      <<" second(s)."<< endl;

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    Info << "Interpolation from cartesian to rotated cartesian mesh...\n" << endl;

    // Iterate for each sun ray
    int iter = 0;

    forAll(sunPosVector, vectorID)
    {
        // start clock
        tstartStep = std::clock();

        ////////////////////////////////////////////////////////////////////////
        // Setup for each sun ray

        // Reference
        vectorList &qrsw = qrswList[vectorID];

        // Initialize LAI
        qrsw = zeroListVectors_nMeshCells;

        // sunPosVector i
        vector n2 = sunPosVector[vectorID];
        n2 /= mag(n2);

        // only if sun is above the horizon
        if (n2[1] > 0)
        {

            //clock_t tstartlocal = std::clock();
            DynamicField<point> startListCells(nMeshCells);
            DynamicField<point> endListCells(nMeshCells);
            List<pointIndexHit> pHitListCells(nMeshCells);

            // Calculate solar short-wave radiation vector field
            forAll(qrsw, cellI)
            {
                point starti = pmeshC[cellI];
                point endi = calcEndPoint(starti, n2, pminO, pmaxO);
                startListCells.append(starti);
                endListCells.append(endi);
            }

            surfacesMesh.findLine(startListCells, endListCells, pHitListCells);

            forAll(qrsw, cellI)
            {
                if (!pHitListCells[cellI].hit())
                {
                    qrsw[cellI] = -n2*IDN[vectorID];
                }
            }

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //iter +=1;
        }

        // Export --- temp.
        if (args.optionFound("writeFields")==true)
        {
            std::ostringstream ss;
            ss << 172800+iter*600;
            Info << "Info: Exporting step " << vectorID << endl;

            volVectorField qrswi
            (
                IOobject
                (
                   "qrsw",
                   word(ss.str()),
                   mesh,
                   IOobject::NO_READ
                ),
                mesh,
                dimensionedVector("0", dimensionSet(1,0,-2,0,0,0,0), vector::zero)//,
                //boundaryTypes
            );

            forAll(qrswi, cellI)
            {
                //LAIi[cellI] = LAI[cellI];
                qrswi[cellI] = qrsw[cellI];
            }
            //LAIi.correctBoundaryConditions();
            qrswi.correctBoundaryConditions();
            qrswi.write();
            //runTime++;
            iter +=1;

            Info << "time = " << word(ss.str()) << endl;
        }

        ////////////////////////////////////////////////////////////////////////
        // Status Info
        Info << "Solar ray direction " << vectorID
             << ", It took " << (std::clock()-tstartStep) / (double)CLOCKS_PER_SEC
             << " second(s)."<< endl;


    }

    // // Status Info
    Info << "\nTotal time took: "
         << (std::clock()-tstart) / (double)CLOCKS_PER_SEC
         <<" second(s).\n"<< endl;

    Info<< "End\n" << endl;
    return 0;

}


// ************************************************************************* //
