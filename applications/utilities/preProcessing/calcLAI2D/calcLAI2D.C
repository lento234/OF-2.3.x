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
    calcLAI by Lento Manickathan, May

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

        /*
        Info<< "    Calculating LAI" << endl;
        volScalarField LAI(IOobject("LAI",
                                    runTime.timeName(),
                                    mesh,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE),
                           mesh,
                           dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0));
       */

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

       Info <<  theta << endl; //.subDict("interpolationSchemes")
       // Define rotation vectors
       vector n1(1,0,0);
       vector n2(cosTheta,0,sinTheta);
       //vector n2(1,0,0);

       // Define rotation matrix
       tensor T(rotationTensor(n2,n1));
       tensor Tinv(rotationTensor(n1,n2));

       // Info << "T: " << T << endl;
       // Info << "Tinv: " << Tinv << endl;


       /////////////// Determine the properties of original coordinate systems

       // Define mesh bounding box
       treeBoundBox allBb(mesh.points());

       // Define mesh centroid
       // point pCentroid = allBb.max()-allBb.min();

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

       // Info << "Info: " << pmin << endl;
       // Info << "Info: " << pmax << endl;


       /////////////// Interpolate from polyMesh to cartesian

       // Define interpolator
       dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes");
       autoPtr<interpolation<scalar> > LAD_interpolator = interpolation<scalar>::New(interpolationDict, LAD);


       /////////////// Define cartesian interpolation grid

       // Cartesian mesh resolution
       point dp(0.01,0,0.01);

       // Increase tolerance of cartesian grid
       pmin.x() -= dp.x();
       pmin.z()  = gMin(pmeshCRot).z();
       //pmin.z() -= 2*dp.z();
       pmax.x() += dp.x();
       pmax.z() += dp.z();

       // Mesh size
       int nx = ceil((pmax.x()-pmin.x())/dp.x());
       int nz = ceil((pmax.z()-pmin.z())/dp.z());

       Info << "Info: " << nx << " " << nz << endl;


       // Generate cartesian grid, LAD and LAI field
       pointField pInterp(nx*nz, point::zero);
       scalarField LADInterp(nx*nz,pTraits<scalar>::zero);
       scalarField LAIInterp(nx*nz,pTraits<scalar>::zero);
       //DynamicList<int> cellIndexList;


       Info << "Info: " << nx << " " << nz << endl;

       int cellIndex;

       /////////////// Interpolate LAD onto interpolation mesh

       for (int k=0; k < nz; k++)
       {
         for (int i=0; i < nx; i++)
         {
           // x,y,z coordinates
           pInterp[k*nx+i].x() = pmin.x() + i*dp.x();
           pInterp[k*nx+i].z() = pmin.z() + k*dp.z();

           // Intersecting cellIndex
           //cellIndex = mesh.findCell(pInterp[k*nx+i]);
           ptemp = transform(Tinv,pInterp[k*nx+i]);

           if ( (ptemp > allBb.min()) && (ptemp < allBb.max()) )
           {
             //cellIndex = mesh.findCell(ptemp); // fast
             //cellIndex = ms.findCell(ptemp,-1,true); slow
             //cellIndex = ms.findCell(ptemp,-1,false); slower
             //cellIndex = ms.findCell(ptemp,0,true); // faster
             cellIndex = ms.findNearestCell(ptemp,0,true); // fastest // most likely handels holes

             if (cellIndex != -1)
             {
               // cellIndexList.append(cellIndex);
               // Interpolate onto scalar field
               // LADInterp[k*nx+i] = LAD[cellIndex]; // nearest neighbour
               LADInterp[k*nx+i] = LAD_interpolator->interpolate(pInterp[k*nx+i],cellIndex); // interpolate
             }

           }

         }
       }

       // Info << "Info: " << LADInterp << endl;

       /////////////// Integrate LAD

       for (int i=0; i < nx; i++)
       {
         for (int k=(nz-2); k>=0; k--)
         {
           LAIInterp[k*nx+i] = LAIInterp[(k+1)*nx+i] + 0.5*(LADInterp[k*nx+i]+LADInterp[(k+1)*nx+i])*dp.z();
         }
       }

       // Info << "Info: " << LAIInterp << endl;


       /////////////// Interpolate LAI from cartesian to original grid

       double xp,zp,xd,zd,c1,c2;
       int i0, k0;

       forAll(LAD, cellI)
       {
            //p = mesh.C()[cellI];
            // Cell center point (rotated coordinate system)
            ptemp = pmeshCRot[cellI];

            //if (LAD[cellI] > 10*SMALL)
            if ( (ptemp.x() >= pmin.x()) && (ptemp.x() <= pmax.x()) && (ptemp.z() >= pmin.z()) && (ptemp.z() <= pmax.z()) )
            {

                // Offset from p min.
                xp = ptemp.x()-pmin.x();
                zp = ptemp.z()-pmin.z();

                // Determine index of lower bound
                i0 = floor(xp/dp.x());
                k0 = floor(zp/dp.z());

                // Lowest point interpolation
                // LAI[cellI] = LAIInterp[k0*nx+i0];

                // Bilinear interpolation

                xd = (ptemp.x()-pInterp[k0*nx+i0].x())/dp.x();
                zd = (ptemp.z()-pInterp[k0*nx+i0].z())/dp.z();

                c1 = LAIInterp[k0*nx+i0]*(1.0-xd) + LAIInterp[k0*nx+(i0+1)]*xd;
                c2 = LAIInterp[(k0+1)*nx+i0]*(1.0-xd) + LAIInterp[(k0+1)*nx+(i0+1)]*xd;

                LAI[cellI] =  c1*(1.0-zd) + c2*zd;


            }
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
