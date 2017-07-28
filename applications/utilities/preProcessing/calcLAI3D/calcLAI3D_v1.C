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
    calcLAI3D by Lento Manickathan, May

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "polyMesh.H"
//#include "meshTools.H"
#include "meshSearch.H"
#include "interpolation.H"

#include <ctime>

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

       if (gSum(LAD) < 10*SMALL)
       {
         Info << "\n\n\nNo vegetation !!\n\n\n" << endl;
       }


       /////////////// tic
       clock_t tstart = std::clock();


       /////////////// Define solar angle rotationTensor

       // Define Solar angles

       scalar PI = 3.14159265359;
       dimensionedScalar theta = vegetationProperties.lookup("phi");//(0-90)*(PI/180);
       scalar cosTheta = cos((theta.value()-90)*(PI/180));
       scalar sinTheta = sin((theta.value()-90)*(PI/180));


       // Define rotation vectors
       vector n1(1,0,0);
       vector n2(cosTheta,0,sinTheta);
       //vector n2(1,1,1);
       //n2 /= mag(n2); // unit vector

       // Define rotation matrix
       tensor T(rotationTensor(n2,n1));
       tensor Tinv(rotationTensor(n1,n2));

       /////////////// Determine the properties of original coordinate systems

       // Define mesh bounding box
       treeBoundBox allBb(mesh.points());

       // Mesh cell centers
       pointField pmeshC = mesh.C();

       // Define search mesh
       meshSearch ms(mesh);

       /////////////// Determine rotated coordinate system mesh
       pointField pmeshCRot = transform(T,pmeshC);

       /////////////// Determine bbox of vegetation (rotated coordinate system)

       point pmin=gMax(pmeshCRot);
       point pmax=gMin(pmeshCRot);
       point ptemp;

       forAll(LAD, cellI)
       {
            if (LAD[cellI] > 10*SMALL)
            {
                ptemp = pmeshCRot[cellI];
                pmin = min(pmin,ptemp);
                pmax = max(pmax,ptemp);
            }
       }

       /////////////// Interpolate from polyMesh to cartesian

       // Define interpolator
       dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes");
       autoPtr<interpolation<scalar> > LAD_interpolator = interpolation<scalar>::New(interpolationDict, LAD);


       /////////////// Define cartesian interpolation grid

       // Cartesian mesh resolution
       point dp(0.2,0.2,0.2);

       // Increase tolerance of cartesian grid
       pmin -= dp;
       pmax += dp;

       // Increase lower bound
       pmin.z()  = gMin(pmeshCRot).z();

       // Mesh size
       int nx = ceil( (pmax.x()-pmin.x()) / dp.x());
       int ny = ceil( (pmax.y()-pmin.y()) / dp.x());
       int nz = ceil( (pmax.z()-pmin.z()) / dp.x());

       // Generate cartesian grid, LAD and LAI field
       pointField pInterp(nx*ny*nz, point::zero);
       scalarField LADInterp(nx*ny*nz,pTraits<scalar>::zero);
       scalarField LAIInterp(nx*ny*nz,pTraits<scalar>::zero);

       Info << "Mesh size: " << nx << "x" << ny << "x" << nz << " : " << pInterp.size() << endl;

       int cellIndex;
       int pIndex;

       /////////////// Interpolate LAD onto interpolation mesh

       for (int k=0; k < nz; k++)
       {
         for (int j=0; j < ny; j++)
         {
           for (int i=0; i < nx; i++)
           {
             // p interp index
             pIndex = (nx*ny)*k + j*nx + i;

             // x,y,z coordinates
             pInterp[pIndex].x() = pmin.x() + i*dp.x();
             pInterp[pIndex].y() = pmin.y() + j*dp.y();
             pInterp[pIndex].z() = pmin.z() + k*dp.z();

             // Intersecting cellIndex
             ptemp = transform(Tinv, pInterp[pIndex]);

             if ( (ptemp > allBb.min()) && (ptemp < allBb.max()) )
             {
               // Find cell
               //cellIndex = mesh.findCell(ptemp); // fast
               //cellIndex = ms.findCell(ptemp,-1,true); // slow
               //cellIndex = ms.findCell(ptemp,-1,false); // slower
               cellIndex = ms.findCell(ptemp,0,true); // faster
               //cellIndex = ms.findNearestCell(ptemp,0,true); // fastest // most likely handels holes
               //Info << "cellIndex: " << cellIndex << " , size: " << mesh.C().size() << endl;

               if (cellIndex != -1)
               {
                 // Interpolate onto scalar field
                 //LADInterp[pIndex] = LAD[cellIndex]; // nearest neighbour
                 LADInterp[pIndex] = LAD_interpolator->interpolate(pInterp[pIndex],cellIndex); // interpolate
                 //Info << "pIndex: " << pIndex << endl;
               }
               else {
                 //Info << "outside, cellIndex: " << cellIndex << endl;
                 LADInterp[pIndex] = -1e99;
               }

             }

           }
         }
       }

       /////////////// Integrate LAD
       int pIndexkp1;
       for (int i=0; i < nx; i++)
       {
         for (int j=0; j < ny; j++)
         {
           for (int k=(nz-2); k>=0; k--)
           {
             pIndex = (nx*ny)*k + j*nx + i; // lower index
             pIndexkp1 = (nx*ny)*(k+1) + j*nx + i; // upper index
             // Finite difference, Euler integration

             if (LADInterp[pIndex] >= 0)
                 LAIInterp[pIndex] = LAIInterp[pIndexkp1] + 0.5*(LADInterp[pIndex]+LADInterp[pIndexkp1])*dp.z();
           }
         }
       }

       /////////////// Interpolate LAI from cartesian to original grid

       int i0, j0, k0;
       double xp,yp,zp;
       double xd,yd,zd;
       double c00, c01, c10, c11, c0, c1;

       forAll(LAD, cellI)
       {
            // Cell center point (rotated coordinate system)
            ptemp = pmeshCRot[cellI];

            //if ( (ptemp.x() >= pmin.x()) && (ptemp.x() <= pmax.x()) && (ptemp.z() >= pmin.z()) && (ptemp.z() <= pmax.z()) )
            //if ((ptemp >= pmin) && (ptemp <= pmax))
            //if (LAD[cellI] > 10*SMALL)
            //if (LAD[cellI] > 10*SMALL)
            //if ((ptemp > pmin) && (ptemp < pmax))
            //if (LAD[cellI] > 10*SMALL)
            if ( (ptemp.x() >= pmin.x()) && (ptemp.x() <= pmax.x()) &&  (ptemp.y() >= pmin.y()) && (ptemp.y() <= pmax.y()) && (ptemp.z() >= pmin.z()) && (ptemp.z() <= pmax.z()) )
            {

                // Offset from p min.
                xp = ptemp.x()-pmin.x();
                yp = ptemp.y()-pmin.y();
                zp = ptemp.z()-pmin.z();

                // Determine index of lower bound
                i0 = floor(xp/dp.x());
                j0 = floor(yp/dp.y());
                k0 = floor(zp/dp.z());

                // Lowest point interpolation
                //LAI[cellI] = LAIInterp[(nx*ny)*k0 + j0*nx + i0];

                // Bilinear interpolation
                xd = (ptemp.x()-pInterp[(nx*ny)*k0 + j0*nx + i0].x())/dp.x();
                yd = (ptemp.y()-pInterp[(nx*ny)*k0 + j0*nx + i0].y())/dp.x();
                zd = (ptemp.z()-pInterp[(nx*ny)*k0 + j0*nx + i0].z())/dp.z();

                // Interpolation in x-dir
                c00 = LAIInterp[(nx*ny)*k0 + j0*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*k0 + j0*nx + i0+1]*xd;
                c01 = LAIInterp[(nx*ny)*(k0+1) + j0*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*(k0+1) + j0*nx + i0+1]*xd;
                c10 = LAIInterp[(nx*ny)*k0 + (j0+1)*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*k0 + (j0+1)*nx + i0+1]*xd;
                c11 = LAIInterp[(nx*ny)*(k0+1) + (j0+1)*nx + i0]*(1-xd) + LAIInterp[(nx*ny)*(k0+1) + (j0+1)*nx + i0+1]*xd;

                // Interpolation in y-dir
                c0 = c00*(1.0-yd) + c10*yd;
                c1 = c01*(1.0-yd) + c11*yd;

                // Interpolate in z-dir
                LAI[cellI] =  c0*(1.0-zd) + c1*zd;

            }
       }

       forAll(LAI, cellI)
       {
         if (LAI[cellI] < 0)
             LAI[cellI] = 0.0;
       }

       LAI.correctBoundaryConditions();


       Info << "It took "<< (std::clock()-tstart) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;

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
        Info<< "    No LAD" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
