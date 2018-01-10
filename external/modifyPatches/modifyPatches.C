/*---------------------------------------------------------------------------*\
   II   II        II  II   Leichtweiss-Institute for Hydraulics          
  II    II  II  II   II    Dep. Hydromechanics and Coastal Eng. 
 II     IIIIIIII    II     Developed by: Hisham El Safti
IIIIII  II  II     II      Email: hsafti@gmail.com
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
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

Description

This utility is modifies "constant/polyMesh/boundary" to change patch types
to empty, wall, cyclic, processor, symmetry or wedge (for use with gmshToFoam)
and to add null patches at end of file (for use with splitMesh) and remove 
existing null patches

Developed by: Hisham El Safti hsafti@gmail.com, January 2013

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "wallPolyPatch.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "primitiveFacePatch.H"
#include "repatchPolyTopoChanger.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


word typeNameWithoutDash (string currentType)
{
    currentType.erase (0,1);
    return word(currentType);
}


void modifyPatchType(Foam::polyMesh& mesh, word name, word patchType)
{
    const label patchI = mesh.boundaryMesh().findPatchID(name);

    if (patchI == -1)
    {
        FatalErrorIn("modifyPatchType(const polyMesh&, const word&, const word&)")
            << "Cannot find patch " << name << " to make it " << patchType << nl
            << "It should be present to be modified" << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    Info << "Modifying patch: " << name << " to type: " << patchType << endl;
   
    repatchPolyTopoChanger repatcher(mesh);
    List<polyPatch*> newPatchPtrList(mesh.boundaryMesh().size());

    forAll(mesh.boundaryMesh(), patchII)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchII];
        if (patchII != patchI)
        {
            if(mesh.boundaryMesh().types()[patchII] == emptyPolyPatch::typeName) 
            {
                const emptyPolyPatch modifiedPatch 
                    (
                        patch.name(),
                        patch.size(),
                        patch.start(),
                        patchII,
                        mesh.boundaryMesh()// polyBoundaryMesh
                    );

                newPatchPtrList[patchII] = modifiedPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();                
            }
            else if(mesh.boundaryMesh().types()[patchII] == wedgePolyPatch::typeName) 
            {
                const wedgePolyPatch modifiedPatch 
                    (
                        patch.name(),
                        patch.size(),
                        patch.start(),
                        patchII,
                        mesh.boundaryMesh()// polyBoundaryMesh
                    );

                newPatchPtrList[patchII] = modifiedPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if(mesh.boundaryMesh().types()[patchII] == cyclicPolyPatch::typeName) 
            {
                const cyclicPolyPatch modifiedPatch 
                    (
                        patch.name(),
                        patch.size(),
                        patch.start(),
                        patchII,
                        mesh.boundaryMesh()// polyBoundaryMesh
                    );

                newPatchPtrList[patchII] = modifiedPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if(mesh.boundaryMesh().types()[patchII] == processorPolyPatch::typeName) 
            {
                const processorPolyPatch modifiedPatch 
                    (
                        patch.name(),
                        patch.size(),
                        patch.start(),
                        patchII,
                        mesh.boundaryMesh(), 0, 1// polyBoundaryMesh
                    );

                newPatchPtrList[patchII] = modifiedPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if(mesh.boundaryMesh().types()[patchII] == symmetryPolyPatch::typeName) 
            {
                const symmetryPolyPatch modifiedPatch 
                    (
                        patch.name(),
                        patch.size(),
                        patch.start(),
                        patchII,
                        mesh.boundaryMesh()// polyBoundaryMesh
                    );

                newPatchPtrList[patchII] = modifiedPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if(mesh.boundaryMesh().types()[patchII] == wallPolyPatch::typeName) 
            {
                const wallPolyPatch modifiedPatch 
                    (
                        patch.name(),
                        patch.size(),
                        patch.start(),
                        patchII,
                        mesh.boundaryMesh()// polyBoundaryMesh
                    );

                newPatchPtrList[patchII] = modifiedPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else
            {
                newPatchPtrList[patchII] = patch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
        }
        else
        {
            if (patchType == wallPolyPatch::typeName)
            {
                const wallPolyPatch wallPatch (
                    name,
                    mesh.boundaryMesh()[patchI].size() ,
                    mesh.boundaryMesh()[patchI].start() ,
                    mesh.boundaryMesh()[patchI].index() , 
                    mesh.boundaryMesh()
                );

                newPatchPtrList[patchII] = wallPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if (patchType == processorPolyPatch::typeName)
            {
                const processorPolyPatch processorPatch (
                    name,
                    mesh.boundaryMesh()[patchI].size() ,
                    mesh.boundaryMesh()[patchI].start() ,
                    mesh.boundaryMesh()[patchI].index() , 
                    mesh.boundaryMesh(), 0, 1
                );

                newPatchPtrList[patchII] = processorPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if (patchType == emptyPolyPatch::typeName)
            {
                const emptyPolyPatch emptyPatch (
                    name,
                    mesh.boundaryMesh()[patchI].size() ,
                    mesh.boundaryMesh()[patchI].start() ,
                    mesh.boundaryMesh()[patchI].index() , 
                    mesh.boundaryMesh()
                );

                newPatchPtrList[patchII] = emptyPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if (patchType == wedgePolyPatch::typeName)
            {
                const wedgePolyPatch wedgePatch (
                    name,
                    mesh.boundaryMesh()[patchI].size() ,
                    mesh.boundaryMesh()[patchI].start() ,
                    mesh.boundaryMesh()[patchI].index() , 
                    mesh.boundaryMesh()
                );

                newPatchPtrList[patchII] = wedgePatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if (patchType == cyclicPolyPatch::typeName)
            {
                const cyclicPolyPatch cyclicPatch (
                    name,
                    mesh.boundaryMesh()[patchI].size() ,
                    mesh.boundaryMesh()[patchI].start() ,
                    mesh.boundaryMesh()[patchI].index() , 
                    mesh.boundaryMesh()
                );

                newPatchPtrList[patchII] = cyclicPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            else if (patchType == symmetryPolyPatch::typeName)
            {
                const symmetryPolyPatch symmetryPatch (
                    name,
                    mesh.boundaryMesh()[patchI].size() ,
                    mesh.boundaryMesh()[patchI].start() ,
                    mesh.boundaryMesh()[patchI].index() , 
                    mesh.boundaryMesh()
                );

                newPatchPtrList[patchII] = symmetryPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchII,
                        patch.size(),
                        patch.start()
                    ).ptr();
            }
            
        }
                
    }
    
    repatcher.changePatches(newPatchPtrList);    
    
    repatcher.repatch();

}






// Add null patch
void addNullPatch(Foam::polyMesh& mesh, const word& name)
{
    const label patchI = mesh.boundaryMesh().findPatchID(name);


    if (patchI != -1 && mesh.boundaryMesh()[patchI].size())
    {
        FatalErrorIn("addNullPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is present but non-zero size"
            << exit(FatalError);
    }

    if (patchI == -1)
    {
      Info << "Adding null patch: " << name << endl;

      repatchPolyTopoChanger repatcher(mesh);
      List<polyPatch*> newPatchPtrList((mesh.boundaryMesh().size() + 1));
      forAll(mesh.boundaryMesh(), patchII)
      {
          const polyPatch& patch = mesh.boundaryMesh()[patchII];

          if(mesh.boundaryMesh().types()[patchII] == emptyPolyPatch::typeName) 
          {
              const emptyPolyPatch modifiedPatch 
                  (
                      patch.name(),
                      patch.size(),
                      patch.start(),
                      patchII,
                      mesh.boundaryMesh()// polyBoundaryMesh
                  );

              newPatchPtrList[patchII] = modifiedPatch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();                
          }
          else if(mesh.boundaryMesh().types()[patchII] == wedgePolyPatch::typeName) 
          {
              const wedgePolyPatch modifiedPatch 
                  (
                      patch.name(),
                      patch.size(),
                      patch.start(),
                      patchII,
                      mesh.boundaryMesh()// polyBoundaryMesh
                  );

              newPatchPtrList[patchII] = modifiedPatch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
          }
          else if(mesh.boundaryMesh().types()[patchII] == cyclicPolyPatch::typeName) 
          {
              const cyclicPolyPatch modifiedPatch 
                  (
                      patch.name(),
                      patch.size(),
                      patch.start(),
                      patchII,
                      mesh.boundaryMesh()// polyBoundaryMesh
                  );

              newPatchPtrList[patchII] = modifiedPatch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
          }
          else if(mesh.boundaryMesh().types()[patchII] == processorPolyPatch::typeName) 
          {
              const processorPolyPatch modifiedPatch 
                  (
                      patch.name(),
                      patch.size(),
                      patch.start(),
                      patchII,
                      mesh.boundaryMesh(), 0, 1// polyBoundaryMesh
                  );

              newPatchPtrList[patchII] = modifiedPatch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
          }
          else if(mesh.boundaryMesh().types()[patchII] == symmetryPolyPatch::typeName) 
          {
              const symmetryPolyPatch modifiedPatch 
                  (
                      patch.name(),
                      patch.size(),
                      patch.start(),
                      patchII,
                      mesh.boundaryMesh()// polyBoundaryMesh
                  );

              newPatchPtrList[patchII] = modifiedPatch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
          }
          else if(mesh.boundaryMesh().types()[patchII] == wallPolyPatch::typeName) 
          {
              const wallPolyPatch modifiedPatch 
                  (
                      patch.name(),
                      patch.size(),
                      patch.start(),
                      patchII,
                      mesh.boundaryMesh()// polyBoundaryMesh
                  );

              newPatchPtrList[patchII] = modifiedPatch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
          }
          else
          {
              newPatchPtrList[patchII] = patch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
          }
      }
      
      const polyPatch patch 
          (
              name,
              0,                  // size
              mesh.nFaces(),      // start
              mesh.boundaryMesh().size(),  // index
              mesh.boundaryMesh()// polyBoundaryMesh
          );
          
      newPatchPtrList[mesh.boundaryMesh().size()] = patch.clone
          (
              mesh.boundaryMesh(),
              mesh.boundaryMesh().size(),
              patch.size(),
              patch.start()
          ).ptr();

      repatcher.changePatches(newPatchPtrList);    

      repatcher.repatch();
    }
}




// Add null patch
void removeNullPatch(Foam::polyMesh& mesh, const word& name)
{
    const label patchI = mesh.boundaryMesh().findPatchID(name);


    if (patchI == -1 || mesh.boundaryMesh()[patchI].size())
    {
        FatalErrorIn("removeNullPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is not present or of non-zero size"
            << exit(FatalError);
    }

    if (patchI != -1)
    {
      Info << "Removing null patch: " << name << endl;

      repatchPolyTopoChanger repatcher(mesh);
      List<polyPatch*> newPatchPtrList((mesh.boundaryMesh().size() - 1));

      label oneLess(0);

      forAll(mesh.boundaryMesh(), patchII)
      {
      
          if (patchI == patchII)
          {
              oneLess=-1;
          } 
          else
          {
              const polyPatch& patch = mesh.boundaryMesh()[patchII];    
              if(mesh.boundaryMesh().types()[patchII] == emptyPolyPatch::typeName) 
              {
                  const emptyPolyPatch modifiedPatch 
                      (
                          patch.name(),
                          patch.size(),
                          patch.start(),
                          patchII+oneLess,
                          mesh.boundaryMesh()// polyBoundaryMesh
                      );

                  newPatchPtrList[patchII+oneLess] = modifiedPatch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();                
              }
              else if(mesh.boundaryMesh().types()[patchII] == wedgePolyPatch::typeName) 
              {
                  const wedgePolyPatch modifiedPatch 
                      (
                          patch.name(),
                          patch.size(),
                          patch.start(),
                          patchII+oneLess,
                          mesh.boundaryMesh()// polyBoundaryMesh
                      );

                  newPatchPtrList[patchII+oneLess] = modifiedPatch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
              }
              else if(mesh.boundaryMesh().types()[patchII] == cyclicPolyPatch::typeName) 
              {
                  const cyclicPolyPatch modifiedPatch 
                      (
                          patch.name(),
                          patch.size(),
                          patch.start(),
                          patchII+oneLess,
                          mesh.boundaryMesh()// polyBoundaryMesh
                      );

                  newPatchPtrList[patchII+oneLess] = modifiedPatch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
              }
              else if(mesh.boundaryMesh().types()[patchII] == processorPolyPatch::typeName) 
              {
                  const processorPolyPatch modifiedPatch 
                      (
                          patch.name(),
                          patch.size(),
                          patch.start(),
                          patchII+oneLess,
                          mesh.boundaryMesh(), 0, 1// polyBoundaryMesh
                      );

                  newPatchPtrList[patchII+oneLess] = modifiedPatch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
              }
              else if(mesh.boundaryMesh().types()[patchII] == symmetryPolyPatch::typeName) 
              {
                  const symmetryPolyPatch modifiedPatch 
                      (
                          patch.name(),
                          patch.size(),
                          patch.start(),
                          patchII+oneLess,
                          mesh.boundaryMesh()// polyBoundaryMesh
                      );

                  newPatchPtrList[patchII+oneLess] = modifiedPatch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
              }
              else if(mesh.boundaryMesh().types()[patchII] == wallPolyPatch::typeName) 
              {
                  const wallPolyPatch modifiedPatch 
                      (
                          patch.name(),
                          patch.size(),
                          patch.start(),
                          patchII+oneLess,
                          mesh.boundaryMesh()// polyBoundaryMesh
                      );

                  newPatchPtrList[patchII+oneLess] = modifiedPatch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
              }
              else
              {
                  newPatchPtrList[patchII+oneLess] = patch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
              }
          }
      }
      
      repatcher.changePatches(newPatchPtrList);    

      repatcher.repatch();


    }
    
}



// Main program:

int main(int argc, char *argv[])
{

    argList::noParallel();
    argList::addOption("empty", "patchName", "change type of patch to empty");
    argList::addOption("wall", "patchName", "change type of patch to wall");
    argList::addOption("cyclic", "patchName", "change type of patch to cyclic");
    argList::addOption("wedge", "patchName", "change type of patch to wedge");
    argList::addOption("symmetry", "patchName", "change type of patch to symmetry");
    argList::addOption("processor", "patchName", "change type of patch to processor");
    argList::addOption("null", "patchName", "creates a new patch of zero size");
    argList::addOption("remove", "patchName", "removes a zero size patch");

    argList::removeOption("doc");
    argList::removeOption("srcDoc");

    argList::addNote("Description: modifyPatches is an application developed by Hisham El Safti");
    argList::addNote("hsafti@gmial.com in January 2013");
    argList::addNote("Purpose: Manipulates the constant/polyMesh/boundary file to modify patches");
    argList::addNote("from patch type to wall, empty, wedge, cyclic, symmetryPlane or processor. ");
    argList::addNote("This is helpful to use with gmshToFoam. The utility allows addition and removal");
    argList::addNote("of null patches for use with splitMesh.\n ");
    argList::addNote("Examples of use:\n");
    argList::addNote("$ modifyPatches -empty frontAndBack \n");
    argList::addNote("$ modifyPatches -null Master \n");
    argList::addNote("Notes:");
    argList::addNote("* The program overwrites constant/polyMesh/boundary ");
    argList::addNote("* Each option takes one argument so it is recommended to ");
    argList::addNote("  use the application once for each patch manipulation to ");
    argList::addNote("  be sure");
    argList::addNote("* A patch modified to processor will always have a 0 processor");
    argList::addNote("  and a 1 neighbour processor");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"

    runTime.setTime(instant(runTime.constant()), 0);
    
    word currentType;

    if (argc == 1)
    {
        Info << "\nPlease provide command line arguments or type "
            << "\"modifyPatches -help\" for help\n" << endl; 
        return 0;
    }

    for (label i = 1; i < argc; i++)
    {
        if (word(argv[i]) == "-empty"         ||
        word(argv[i]) == "-wedge"             ||
        word(argv[i]) == "-cyclic"            ||
        word(argv[i]) == "-wall"              ||
        word(argv[i]) == "-symmetryPlane"     ||  
        word(argv[i]) == "-processor"         ||
        word(argv[i]) == "-null"              ||
        word(argv[i]) == "-remove"
        )  
        {
            currentType = typeNameWithoutDash((word(argv[i])).c_str());
        }
        else
        {
            if (currentType == "null")
            {
                addNullPatch(mesh, word(argv[i]));
            }
            else if (currentType == "remove")
            {
                removeNullPatch(mesh, word(argv[i]));
            }
            else
            {
                modifyPatchType(mesh, word(argv[i]), currentType);
            }
                
            if (!mesh.write())
            {
                FatalErrorIn(args.executable()) << "Failed writing mesh"
                    << exit(FatalError);
            }
        } 
    }


    Info << "\nModification of patches completed successfully!\n" << endl;     
        
    return 0;

}

// ************************************************************************* //
