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
    initLAD

Description
    initialize LAD by Lento Manickathan.


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

scalar nearestNeighbor3D
(
    const point &iptarget,
    const pointField &psource,
    const scalarIOList &asource,
    const label &nx,
    const label &ny
)
{
    point pmin = psource[0];
    point dp = psource[(nx*ny) + nx + 1] -psource[0]; 

    point p = iptarget - pmin;
    
    // Determine index of lower bound
    int i0 = floor(p.x()/dp.x());
    int j0 = floor(p.y()/dp.y());
    int k0 = floor(p.z()/dp.z()); // if (k0 < 0) k0 = 0;
    
    int i000 = (nx*ny)*k0 + j0*nx + i0;
    Info << "i0: " << i000 << endl;

    return asource[i000];
}


// trilinear interpolation
scalar interp3D
(
    const point &ptemp,
    const pointField &pInterp,
    const scalarIOList &valInterp,
    const label &nx,
    const label &ny
)
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
  int k0 = floor(zp/dp.z()); // if (k0 < 0) k0 = 0;

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    // start timer
    clock_t tstart = std::clock();

    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    Info << "timeDirs: " << timeDirs << endl;
    runTime.setTime(timeDirs[0], 0);

    #include "createNamedMesh.H"

    // Target leaf area density field
    volScalarField a 
    (
      IOobject
      (
          "a",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      ),
      mesh
    );
  
    // Load vegetation properties
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
   
    // Cartesian grid properties
    point pmin = vegetationProperties.lookup("pmin");
    point np  = vegetationProperties.lookup("np");
    point dp = vegetationProperties.lookup("dp");
    
    Info << "Cartesian grid properties ..." << endl;
    Info << "pmin = " << pmin << endl;
    Info << "np = " << np << endl;
    Info << "dp = " << dp << endl;
    
    label N = np[0] * np[1] * np[2];

    Info << "N = " << N << endl;
   
    // Load source leaf area density
    scalarIOList a_source 
    (
        IOobject
        (
            "asource",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (a_source.size() != N)
        FatalErrorIn("a_source.size() != N") << "size does not match input: "
            << a_source.size() << " =! " << N
            << abort(FatalError);
    
    ////////////////////////////////////////////////////////////////////////////
    // Global setup
    
    Info << "Initializing leaf area density (a) field ... \n" << endl;

    // Mesh of source
    pointField p_source(N);

    label pindex;
    for (label k=0; k < np[2]; ++k)
    {
      for (label j=0; j < np[1]; ++j)
      {
        for (label i=0; i < np[0]; ++i)
        {
          // index of node p
          pindex = (np[0]*np[1])*k + j*np[0] + i;

          // x,y,z coordinates
          p_source[pindex].x() = pmin.x() + i*dp.x();
          p_source[pindex].y() = pmin.y() + j*dp.y();
          p_source[pindex].z() = pmin.z() + k*dp.z();
        }
      }
    }

    // Bounding box max
    point pmax = gMax(p_source);

    Info << "bounding box size = " << pmax - pmin << endl;

    // Mesh cell centers
    pointField p_target = mesh.C();
    scalar leafarea;
    
    // Loop through target field
    forAll(a, cellI)
    {
        // Point 
        point ip_target = p_target[cellI];

        // Check if inside bounding box
        if ((ip_target.x() > pmin.x() && ip_target.x() < pmax.x()) &&
            (ip_target.y() > pmin.y() && ip_target.y() < pmax.y()) &&
            (ip_target.z() > pmin.z() && ip_target.z() < pmax.z()))
        {
            // Tri-linear interpolation
            //a[cellI] = interp3D(ip_target, p_source, a_source, np[0], np[1]);
            leafarea = interp3D(ip_target, p_source, a_source, np[0], np[1]); //
            //a[cellI] = leafarea/mesh.V()[cellI];
            a[cellI] = leafarea;///mesh.V()[cellI];
            //Info << "Mesh volume: " << meshV()[cellI] << endl;
            //a[cellI] = nearestNeighbor3D(ip_target, p_source, a_source, np[0], np[1]);
        }
    }

    Info << "gMin = " << gMin(a) << endl;
    Info << "gMax = " << gMax(a) << endl;

    // Write out leaf area density field
    //a.correctBoundaryConditions();
    a.write();

    // // Status Info
    Info << "\nTotal time took: "
         << (std::clock()-tstart) / static_cast<double>(CLOCKS_PER_SEC)
         <<" second(s).\n"<< endl;

    Info<< "End\n" << endl;
    return 0;

}


// ************************************************************************* //
