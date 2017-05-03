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
    raytracing

Description
    raytracing by Lento Manickathan, April.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "polyMesh.H"
#include "meshTools.H"
#include "meshSearch.H"
#include "treeDataFace.H"
#include "treeDataCell.H"
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

    if (LADheader.headerOk())
    {
        Info<< "    Reading LAD" << endl;
        volScalarField LAD(LADheader, mesh);

        Info<< "    Calculating LAI" << endl;
        volScalarField LAI(IOobject("LAI",
                                    runTime.timeName(),
                                    mesh,
                                    IOobject::NO_READ,
                                    IOobject::AUTO_WRITE),
                           mesh,
                           dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0));


       // Determine the mesh bounding box
       treeBoundBox allBb(mesh.points());
       //scalar bbTol = 1e-6 * allBb.avgDim();
       //point& bbMin = allBb.min();
       //bbMin.x() -= bbTol;
       //bbMin.y() -= bbTol;
       //bbMin.z() -= bbTol;



       clock_t tstart = std::clock();

       /*
       // Generate octree of mesh cells
       indexedOctree<treeDataCell> cellTree
       (
           treeDataCell(false, mesh, polyMesh::FACECENTRETETS),
           allBb, // overall search domain
           8, // maxLevel
           10, // leafsize
           3.0 // duplicity
       );
       */


       double xmin=1e10;
       double xmax=SMALL;
       double zmin=1e10;
       double zmax=SMALL;
       point p;

       forAll(LAD, cellI)
       {
            if (LAD[cellI] > 10*SMALL)
            {
                p = mesh.C()[cellI];
                xmin = min(xmin,p.x());
                xmax = max(xmax,p.x());
                zmin = min(zmin,p.z());
                zmax = max(zmax,p.z());
            }
       }

       /*
       Info << "xmin: " << xmin << endl;
       Info << "zmin: " << zmin << endl;
       Info << "xmax: " << xmax << endl;
       Info << "zmax: " << zmax << endl;
       */

       // vector<point> v1(Foam::vector::zero);


       dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes");
       autoPtr<interpolation<scalar> > LAD_interpolator = interpolation<scalar>::New(interpolationDict, LAD);

       double dx = 0.01;
       double dz = 0.01;

       // Increase tolerance of bounding domain
       xmin -= 2*dx;
       //zmin -= 2*dz;
       zmin = allBb.min().z();
       xmax += 2*dx;
       zmax += 2*dz;

       int nx = ceil((xmax-xmin)/dx);
       int nz = ceil((zmax-zmin)/dz);

       pointField pInterp(nx*nz, point::zero);
       scalarField LADInterp(nx*nz,pTraits<scalar>::zero);
       scalarField LAIInterp(nx*nz,pTraits<scalar>::zero);
       DynamicList<int> cellIndexList;

       //DynamicList<scalar> pData;
       int cellIndex;

       // Interpolate LAD onto interpolation mesh
       for (int k=0; k < nz; k++)
       {
         for (int i=0; i < nx; i++)
         {
           // x,y,z coordinates
           pInterp[k*nx+i].x() = xmin + i*dx;
           pInterp[k*nx+i].z() = zmin + k*dz;

           // Intersecting cellIndex
           cellIndex = mesh.findCell(pInterp[k*nx+i]);

           if (cellIndex != -1)
           {
             cellIndexList.append(cellIndex);
             // Interpolate onto scalar field
             // LADInterp[k*nx+i] = LAD[cellIndex]; // nearest neighbour
             LADInterp[k*nx+i] = LAD_interpolator->interpolate(pInterp[k*nx+i],cellIndex); // interpolate
           }

         }
       }

       // Integrate LAD
       for (int i=0; i < nx; i++)
       {
         for (int k=(nz-2); k>=0; k--)
         {
           LAIInterp[k*nx+i] = LAIInterp[(k+1)*nx+i] + 0.5*(LADInterp[k*nx+i]+LADInterp[(k+1)*nx+i])*dz;
         }
       }
       //
       //p1[0] = mesh.C()[1000];

       //Info << "Info: " << p1[0] << endl;
       //Info << "Info: " << LADInterp << endl;
       //Info << "Info: " << LAIInterp << endl;

       //Info << "Info: " << cellIndexList << endl;


       double xp,zp;
       int i0, k0;


       forAll(LAD, cellI)
       {
            p = mesh.C()[cellI];

            //if (LAD[cellI] > 10*SMALL)
            if ( (p.x() >= xmin) && (p.x() <= xmax) && (p.z() >= zmin) && (p.z() <= zmax) )
            {
                // Interpolate
                //p = mesh.C()[cellI];

                // Offset from p min.
                xp = p.x()-xmin;
                zp = p.z()-zmin;

                // Determine index of lower bound
                i0 = floor(xp/dx);
                k0 = floor(zp/dz);

                // Lowest point interpolation
                //LAI[cellI] = LAIInterp[k0*nx+i0];

                // Bilinear interpolation

                double xd=(p.x()-pInterp[k0*nx+i0].x())/dx;
                double zd=(p.z()-pInterp[k0*nx+i0].z())/dz;

                double c1=LAIInterp[k0*nx+i0]*(1.0-xd) + LAIInterp[k0*nx+(i0+1)]*xd;
                double c2=LAIInterp[(k0+1)*nx+i0]*(1.0-xd) + LAIInterp[(k0+1)*nx+(i0+1)]*xd;

                LAI[cellI] =  c1*(1.0-zd) + c2*zd;

            }
       }

       Info << "It took "<< (std::clock()-tstart) / (double)CLOCKS_PER_SEC <<" second(s)."<< endl;



       /*
       int k = 0;

       dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes");
       autoPtr<interpolation<scalar> > LADint =interpolation<scalar>::New(interpolationDict, LAD);

       double PI = 3.14159265359;
       double theta = 45*(PI/180);
       double cosTheta = cos(theta);
       double sinTheta = sin(theta);


       pointIndexHit pHit;
       scalar dL;
       scalar interpValue;

       */

          //  forAll(LAD, cellI)
          //  {
          //    // Integral LAI value
          //    if ((LAD[cellI] > 10*SMALL)) // && (LAI[cellI] < SMALL))
          //    {
           //
          //      k++;
           //
          //      // Determine cell center
          //      point pCell = mesh.C()[cellI]; // start point of line
           //
          //      // Solar angle
           //
          //      // double tanTheta = tan(theta);
           //
          //      // Determine starting and ending point // z --> vertical
          //         point pStart(pCell.x(), pCell.y(), pCell.z());
           //
          //      //point pEnd(pCell.x() - 50*cosTheta, pCell.y(), pCell.z() - 50*sinTheta);
           //
          //         point pEnd(pCell.x() + 50*cosTheta, pCell.y(), pCell.z() + 50*sinTheta);
           //
          //      //point pStart(pCell.x() + (allBb.max().z()-pCell.z())/tanTheta, pCell.y(), allBb.max().z());
          //      //point pEnd(pCell.x() - (allBb.min().z()-pCell.z())/tanTheta, pCell.y(), allBb.min().z());
           //
          //      //point pEnd(pCell.x() - 50*cosTheta, pCell.y(), pCell.z() - 50*sinTheta);
          //      //point pEnd(pCell.x() - (allBb.min().z()-pCell.z())/tanTheta, pCell.y(), allBb.min().z());
          //      //point pStart(pCell.x() + (allBb.max().z()-pCell.z())/tanTheta, pCell.y(), allBb.max().z()); // end point of line
           //
           //
          //      // Determine starting and ending point // y --> vertical
          //      //////////////    point pEnd(pCell.x(), pCell.y(), pCell.z());
          //      //point pEnd(pCell.x() - 50*cosTheta, pCell.y(), pCell.z() - 50*sinTheta);
          //      //////////////    point pStart(pCell.x() + 50*cosTheta, pCell.y()  + 50*sinTheta, pCell.z());
          //      //point pStart(pCell.x() + (allBb.max().z()-pCell.z())/tanTheta, pCell.y(), allBb.max().z());
          //      //point pEnd(pCell.x() - (allBb.min().z()-pCell.z())/tanTheta, pCell.y(), allBb.min().z());
           //
          //      //point pEnd(pCell.x() - 50*cosTheta, pCell.y(), pCell.z() - 50*sinTheta);
          //      //point pEnd(pCell.x() - (allBb.min().z()-pCell.z())/tanTheta, pCell.y(), allBb.min().z());
          //      //point pStart(pCell.x() + (allBb.max().z()-pCell.z())/tanTheta, pCell.y(), allBb.max().z()); // end point of line
           //
           //
          //      // Define the direction
          //      const vector eVec(pEnd - pStart); // line vector
          //      const vector tolVec = 1e-10*eVec;//1e-6*eVec;
           //
          //      // Integral LAI value
          //      scalar value = 0.0;
           //
          //      // ray-tracing loop
          //      while (true)
          //      {
          //        // Hit point
          //        //pointIndexHit pHit = cellTree.findLine(pStart, pEnd);
          //        pHit = cellTree.findLine(pStart, pEnd);
           //
          //        // If point is hit
          //        if (pHit.hit())
          //        {
           //
          //         // Set LAI of present cell
          //         //LAI[pHit.index()] = value;
          //         //Info << "Hit point " << pHit.hitPoint() << endl;
           //
          //         // Determine dL, the integrant
          //         //scalar dL = mag(pHit.hitPoint()-pStart);
          //         dL = mag(pHit.hitPoint()-pStart);
           //
          //         // Integrate leaf area density,
          //         // scalar interpValue = LADint->interpolate(pHit.hitPoint(), pHit.index());
          //         interpValue = LADint->interpolate(pHit.hitPoint(), pHit.index());
          //         //Info << "interpValue: " << interpValue << endl;
          //         //Info << "cell center value: " << LAD[pHit.index()] << endl;
          //         //Info << "error: " << mag(interpValue-LAD[pHit.index()]) << endl;
          //         //value += LAD[pHit.index()]*dL; // cell-center value
          //         value += interpValue*dL; // interpolated
           //
          //         /*
          //         // LAI already calculated in the hit point
          //         if (LAI[pHit.index()] > 10*SMALL)
          //         {
          //           break;
          //         }
          //         */
          //         if (LAD[pHit.index()] < 10*SMALL)
          //         {
          //           break;
          //         }
           //
          //         // set new start point shortly after previous start point
          //         pStart = pHit.hitPoint() + tolVec;
           //
          //       }
          //       else
          //       {
          //         // No hit.
          //         // Info << "No hit" << endl;
          //         break;
          //       }
           //
          //     } // end of ray tracing
          //     LAI[cellI] = value;
           //
          //     // Set LAI of index of the last cell
          //     //LAI[cellI] = value;
          //     // Info
          //     // Info << "LAI: " << LAI[cellI] << endl;
          //     Info << "k: " << k << endl;
          //    } // if inside vegetation
          //    // if (k > 100) break;
           //
          //  } // iterate through all cells



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
