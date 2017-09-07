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
    calcLAI3D

Description
    calcLAI3D by Lento Manickathan

Versions
    May   - v1
    July  - v2
    Aug   - v3

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "wallPolyPatch.H"
#include "treeDataFace.H"
#include "vectorIOList.H"
#include <ctime>
//#include "meshTools.H"
//#include "treeDataCell.H"
//#include "treeDataCell.H"


namespace Foam
{

// calculate the end point for a ray hit check
point calcEndPoint(point &start, point &n2, point &pminO, point &pmaxO)
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

/*
// trilinear interpolation
scalar interp3D(scalarField &LAIInterp, pointField &pInterp, point &ptemp,  point &dp, int &i0, int &j0, int &k0, int &nx, int &ny, int &nz)
{

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
  scalar xd = (ptemp.x()-pInterp000.x())/dp.x();
  scalar yd = (ptemp.y()-pInterp000.y())/dp.y();
  scalar zd = (ptemp.z()-pInterp000.z())/dp.z();
  // xd = (ptemp.x()-pInterp[(nx*ny)*k0 + j0*nx + i0].x())/dp.x();
  // yd = (ptemp.y()-pInterp[(nx*ny)*k0 + j0*nx + i0].y())/dp.y();
  // zd = (ptemp.z()-pInterp[(nx*ny)*k0 + j0*nx + i0].z())/dp.z();

  // Interpolation in x-dir
  scalar c00 = LAIInterp[i000]*(1-xd) + LAIInterp[i100]*xd;
  scalar c01 = LAIInterp[i001]*(1-xd) + LAIInterp[i101]*xd;
  scalar c10 = LAIInterp[i010]*(1-xd) + LAIInterp[i110]*xd;
  scalar c11 = LAIInterp[i011]*(1-xd) + LAIInterp[i111]*xd;
  // c00 = LAIInterp[(nx*ny)*k0 + j0*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*k0 + j0*nx + i0+1]*xd;
  // c01 = LAIInterp[(nx*ny)*(k0+1) + j0*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*(k0+1) + j0*nx + i0+1]*xd;
  // c10 = LAIInterp[(nx*ny)*k0 + (j0+1)*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*k0 + (j0+1)*nx + i0+1]*xd;
  // c11 = LAIInterp[(nx*ny)*(k0+1) + (j0+1)*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*(k0+1) + (j0+1)*nx + i0+1]*xd;

  // Interpolation in y-dir
  scalar c0 = c00*(1.0-yd) + c10*yd;
  scalar c1 = c01*(1.0-yd) + c11*yd;

  // Interpolate in z-dir
  return c0*(1.0-zd) + c1*zd;

}
*/

// trilinear interpolation
scalar interp3D(point &ptemp, pointField &pInterp, scalarField &LAIInterp, point &dp, int &nx, int &ny)
{

  // point &dp, int &i0, int &j0, int &k0, int &nx, int &ny, int &nz)
  point pmin = pInterp[0];

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
  double c00 = LAIInterp[i000]*(1-xd) + LAIInterp[i100]*xd;
  double c01 = LAIInterp[i001]*(1-xd) + LAIInterp[i101]*xd;
  double c10 = LAIInterp[i010]*(1-xd) + LAIInterp[i110]*xd;
  double c11 = LAIInterp[i011]*(1-xd) + LAIInterp[i111]*xd;

  // Interpolation in y-dir
  double c0 = c00*(1.0-yd) + c10*yd;
  double c1 = c01*(1.0-yd) + c11*yd;

  // Interpolate in z-dir
  return c0*(1.0-zd) + c1*zd;

}


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
  {
      bool writeResults = !args.optionFound("noWrite");

      IOobject LADheader
      (
          "LAD",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ
      );

      IOobject LAIheader
      (
          "LAI",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ
      );


      if (LADheader.headerOk() && LAIheader.headerOk())
      {
          Info<< "    Reading LAD" << endl;
          volScalarField LAD(LADheader, mesh);

          Info<< "    Reading LAI" << endl;
          volScalarField LAI(LAIheader, mesh);


          IOdictionary vegetationProperties
          (
             IOobject
             (
                 "vegetationProperties",
                 runTime.constant(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::NO_WRITE
             )
          );

          // Check if vegetation is present
          if (gSum(LAD) < 10*SMALL)
          {
           Info << "\n\n\nNo vegetation !!!!!!!!!\n\n\n" << endl;
          }


          //forAll(sunPosVector, vectorId)

          // start timer
          clock_t tstart = std::clock();
          clock_t tstartlocal;

          /////////////// (Step 1.a) Define solar angle
          tstartlocal = std::clock();

          // Define Solar angles (z is vertical)
          //dimensionedScalar phi = vegetationProperties.lookup("phi");    // solar azimuth (degrees)
          //dimensionedScalar beta = vegetationProperties.lookup("beta");  // solar altitude (degrees)

          //dimensionedScalar n2x = vegetationProperties.lookup("n2x");
          //dimensionedScalar n2y = vegetationProperties.lookup("n2y");
          //dimensionedScalar n2z = vegetationProperties.lookup("n2z");
          vector n2 = vegetationProperties.lookup("n2");

          // calculate sines and cosines
          /*
          scalar PI = 3.14159265359;
          scalar cosPhi = cos(phi.value()*(PI/180));
          scalar sinPhi = sin(phi.value()*(PI/180));
          scalar cosBeta = cos(beta.value()*(PI/180));
          scalar sinBeta = sin(beta.value()*(PI/180));
          */

          // Info << "Info: phi = " << phi << ", beta = " << beta << endl;


          /////////////// (Step 1.b) Define rotation matrix

          // Define rotation vectors
          vector n1(0,0,1);                                   // original vector

          //vector n2(cosBeta*cosPhi, cosBeta*sinPhi, sinBeta); // solar vector
          // convert to unit vectors
          n1 /= mag(n1);
          n2 /= mag(n2);

          Info << "Info: n2 = " << n2 << endl;

          // Define rotation matrix
          tensor T(rotationTensor(n2,n1));       // from n1 to n2
          tensor Tinv(rotationTensor(n1,n2));    // from n2 back to n1


          /////////////// (Step 1.c) Determine the properties of mesh

          // Define mesh bounding box
          treeBoundBox allBb(mesh.points());

          // Define search mesh
          meshSearch ms(mesh);

          // Mesh cell centers
          pointField pmeshC = mesh.C();

          // mesh bounding box
          point pminO = gMin(pmeshC);
          point pmaxO = gMax(pmeshC);

          // Define boundary mesh
          const polyBoundaryMesh& patches = mesh.boundaryMesh();

          // Number of boundary faces
          labelList bndFaces(mesh.nFaces()-mesh.nInternalFaces());

          // no idea what this does - magic
          label bndI = 0;
          forAll(patches, patchI)
          {
              const polyPatch& pp = patches[patchI];

              if (!pp.coupled())
              {
                  forAll(pp, i)
                  {
                      bndFaces[bndI++] = pp.start() + i;
                  }
              }
          }
          bndFaces.setSize(bndI);

          // Define boundary mesh octree
          indexedOctree<treeDataFace> boundaryTree
          (
              treeDataFace    // all information needed to search faces
              (
                  false,      // do not cache bb
                  mesh,
                  bndFaces    // boundary faces only
              ),
              allBb,          // overall search domain
              8,              // maxLevel
              10,             // leafsize
              3.0             // duplicity
          );


          /////////////// (Step 1.d) Determine bbox of vegetation (rotated coordinate system)
          // Calculated in the rotated coordinate system

          // Mesh cell centers (rotated coordinate system)
          pointField pmeshCRot = transform(T,pmeshC);

          // Minimum point
          point pmeshMinRot = gMin(pmeshCRot);

          point pmin = gMax(pmeshCRot);
          point pmax = gMin(pmeshCRot);

          int vegetationCell = 0;

          point ptemp;
          forAll(LAD, cellI)
          {
              // where vegetation is present
              if (LAD[cellI] > 10*SMALL)
              {
                  ptemp = pmeshCRot[cellI];
                  pmin = min(pmin,ptemp);
                  pmax = max(pmax,ptemp);
                  vegetationCell = cellI;
              }
          }

          /////////////// (Step 1.e) Define LAD interpolator to arbitrary locations

          // Read interpolation scheme
          dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes");

          // Define interpolator
          autoPtr<interpolation<scalar> > LAD_interpolator = interpolation<scalar>::New(interpolationDict, LAD);

          //////////////////////////////////////////////////////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////////////////


          /////////////// (Step 2) Define cartesian interpolation grid

          // Cartesian mesh resolution (determine from minimum cell size)
          scalar minCellV = gMin(mesh.V());
          scalar minCellL = pow(minCellV, 1.0/3.0);

          // grid spacing
          point dp(minCellL,minCellL,minCellL); //point dp(0.5, 0.5, 0.5);

          // Extend the cartesian grid to include vegetation
          pmin -= 5*dp;
          pmax += 5*dp;


          // Define cartesian grid size
          int nx = ceil( (pmax.x()-pmin.x()) / dp.x());
          int ny = ceil( (pmax.y()-pmin.y()) / dp.y());
          int nz = ceil( (pmax.z()-pmin.z()) / dp.z());

          // Generate cartesian interpolation grid
          // coordinates
          pointField pInterp(nx*ny*nz, point::zero);
          // interpolated LAD
          scalarField LADInterp(nx*ny*nz,pTraits<scalar>::zero);
          // interpolated LAI
          scalarField LAIInterp(nx*ny*nz,pTraits<scalar>::zero);

          // Grid info
          Info << "Info: Mesh size: " << nx << "x" << ny << "x" << nz << " : " << pInterp.size() << endl;

          Info << "Info: <Setup of mesh> took "<< (std::clock()-tstartlocal) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;

          //////////////////////////////////////////////////////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////////////////

          /////////////// (Step 3) Interpolate LAD onto cartesian interpolation mesh
          tstartlocal = std::clock();

          int cellIndex;
          int pIndex;

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
                ptemp = transform(Tinv, pInterp[pIndex]);

                // check if point is inside bounding box of all domain
                //if ( (ptemp > allBb.min()) && (ptemp < allBb.max()) )
                if ( (ptemp >= pmin) && (ptemp <= pmax) )
                {
                  // Find intersecting cell

                  //cellIndex = mesh.findCell(ptemp); // fast
                  //cellIndex = ms.findCell(ptemp,-1,true); // slow
                  //cellIndex = ms.findCell(ptemp,-1,false); // slower
                  //cellIndex = ms.findCell(ptemp,0,true); // faster
                  cellIndex = ms.findCell(ptemp,vegetationCell,true); // faster
                  //cellIndex = ms.findNearestCell(ptemp,0,true); // fastest // most likely handels holes

                  // if point is inside domain
                  if (cellIndex != -1)
                    LADInterp[pIndex] = LAD_interpolator->interpolate(pInterp[pIndex],cellIndex); // interpolate //LADInterp[pIndex] = LAD[cellIndex]; // nearest neighbour

                }

              }
            }
          }

          Info << "Info: <Interpolate LAD to cartesian> took "<< (std::clock()-tstartlocal) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;


          //////////////////////////////////////////////////////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////////////////

          /////////////// (Step 5) Calculated LAI by integrating LAD vertical (from top to bottom)
          tstartlocal = std::clock();

          int pIndexkp1;

          for (int i=0; i < nx; i++)
          {
            for (int j=0; j < ny; j++)
            {
              for (int k=(nz-2); k>=0; k--)
              {
                // lower and upper row index
                pIndex = (nx*ny)*k + j*nx + i;
                pIndexkp1 = (nx*ny)*(k+1) + j*nx + i;

                // trapezoidal integration
                LAIInterp[pIndex] = LAIInterp[pIndexkp1] + 0.5*(LADInterp[pIndex]+LADInterp[pIndexkp1])*dp.z();

              }
            }
          }
          // Info << "gMax(LAIInterp): " << gMax(LAIInterp) << endl;

          Info << "Info: <LAD integration> took "<< (std::clock()-tstartlocal) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;


          //////////////////////////////////////////////////////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////////////////

          /////////////// (Step 6) Interpolate LAI from cartesian to original grid
          tstartlocal = std::clock();

          //int i0, j0, k0;
          //double xp,yp,zp;
          //double xd,yd,zd;
          //double c00, c01, c10, c11, c0, c1;
          point start, end;
          pointIndexHit pHit;


          forAll(LAD, cellI)
          {
              // Cell center point (rotated coordinate system)
              ptemp = pmeshCRot[cellI];

              if ( (ptemp.x() >= pmin.x()) && (ptemp.x() <= pmax.x()) &&
                   (ptemp.y() >= pmin.y()) && (ptemp.y() <= pmax.y()) &&
                   (ptemp.z() >= pmeshMinRot.z()) && (ptemp.z() <= pmax.z()) )
              {

                  // Offset from p min.
                  //xp = ptemp.x()-pmin.x();
                  //yp = ptemp.y()-pmin.y();
                  //zp = ptemp.z()-pmin.z();

                  // Determine index of lower bound
                  //i0 = floor(xp/dp.x());
                  //j0 = floor(yp/dp.y());
                  //k0 = max(0, floor(zp/dp.z())); // if (k0 < 0) k0 = 0;

                  // Check if cell in building shadow shadow
                  start = transform(Tinv, ptemp);
                  end = calcEndPoint(start, n2, pminO, pmaxO);
                  pHit = boundaryTree.findLine(start, end);

                  // If not in the shadow, then set LAI
                  if (pHit.hit() == false)
                  {
                    //LAI[cellI] = interp3D(LAIInterp, pInterp, ptemp, dp, i0, j0, k0, nx, ny, nz);
                    LAI[cellI] = interp3D(ptemp, pInterp, LAIInterp, dp, nx, ny);
                  }
                  else if (LAD[cellI] > 10*SMALL)
                  {
                    LAI[cellI] = 1e3;
                  }

              }
          }

          LAI.correctBoundaryConditions();


          Info << "Info: <LAI interpolation to mesh> took "<< (std::clock()-tstartlocal) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;

          Info << "Info: << Total time took >> "<< (std::clock()-tstart) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;

          LAI.write();


          if (writeResults)
          {
              LAI.write();
          }
          else
          {
              Info<< "        No write"  << endl;
          }

      }
      else
      {
          Info<< "    No LAD or LAI" << endl;
      }

      Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
