/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "directAndDiffuse.H"
#include "surfaceFields.H"
#include "constants.H"
#include "solarLoadViewFactorFixedValueFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "typeInfo.H"
#include "Time.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solarLoad
    {
        defineTypeNameAndDebug(directAndDiffuse, 0);
        addToSolarLoadRunTimeSelectionTables(directAndDiffuse);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarLoad::directAndDiffuse::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::GeometricBoundaryField& Qsp = Qs_.boundaryField();

    label count = 0;
    forAll(Qsp, patchI)
    {
        //const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        const fvPatchScalarField& QsPatchI = Qsp[patchI];

        if ((isA<fixedValueFvPatchScalarField>(QsPatchI)))
        {
            selectedPatches_[count] = QsPatchI.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchI].size();
			
			if ((isA<wallFvPatch>(mesh_.boundary()[patchI])))
			{
				wallPatchOrNot_[count] = 1;
				nLocalWallCoarseFaces_ += coarsePatches[patchI].size();
			}	
			
			count++;
        }
    }
	Info<< "Selected patches:" << selectedPatches_ << endl;
	Info<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
	Info << "wallPatchOrNot_: " << wallPatchOrNot_ << endl;
	Info << "nLocalWallCoarseFaces_: " << nLocalWallCoarseFaces_ << endl;
	
    selectedPatches_.resize(count--);
	wallPatchOrNot_.resize(count--);

	Info<< "Selected patches:" << selectedPatches_ << endl;
	Info<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
	Info << "wallPatchOrNot_: " << wallPatchOrNot_ << endl;
	Info << "nLocalWallCoarseFaces_: " << nLocalWallCoarseFaces_ << endl;

	if (debug)
    {
        Pout<< "Selected patches:" << selectedPatches_ << endl;
        Pout<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());

    if (Pstream::master())
    {
        Info<< "Total number of clusters : " << totalNCoarseFaces_ << endl;
    }

    labelListIOList subMap
    (
        IOobject
        (
            "subMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList constructMap
    (
        IOobject
        (
            "constructMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    IOList<label> consMapDim
    (
        IOobject
        (
            "constructMapDim",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );    

    map_.reset
    (
        new mapDistribute
        (
            consMapDim[0],
            Xfer<labelListList>(subMap),
            Xfer<labelListList>(constructMap)
        )
    );

    scalarListIOList FmyProc
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
	
    scalarListIOList skyViewCoeffmyProc
    (
        IOobject
        (
            "skyViewCoeff",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    skyViewCoeffSize = skyViewCoeffmyProc.size();	
	
    scalarListIOList sunViewCoeffmyProc
    (
        IOobject
        (
            "sunViewCoeff",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    sunViewCoeffSize = sunViewCoeffmyProc.size();

    labelIOList sunskyMapmyProc
    (
        IOobject
        (
            "sunskyMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );    

    labelListIOList globalFaceFaces
    (
        IOobject
        (
            "globalFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ); 

    List<labelList> sunskyMap(Pstream::nProcs());
    sunskyMap[Pstream::myProcNo()] = sunskyMapmyProc;
    Pstream::gatherList(sunskyMap);

    List<labelListList> globalFaceFacesProc(Pstream::nProcs());
    globalFaceFacesProc[Pstream::myProcNo()] = globalFaceFaces;
    Pstream::gatherList(globalFaceFacesProc);

    List<scalarListList> F(Pstream::nProcs());
    F[Pstream::myProcNo()] = FmyProc;
    Pstream::gatherList(F);	
	
    List<scalarListList> skyViewCoeff(Pstream::nProcs());
    skyViewCoeff[Pstream::myProcNo()] = skyViewCoeffmyProc;
    Pstream::gatherList(skyViewCoeff);	
	
    List<scalarListList> sunViewCoeff(Pstream::nProcs());
    sunViewCoeff[Pstream::myProcNo()] = sunViewCoeffmyProc;
    Pstream::gatherList(sunViewCoeff);		

    globalIndex globalNumbering(nLocalCoarseFaces_);

    if (Pstream::master())
    {
        Fmatrix_.reset
        (
            new scalarSquareMatrix(totalNCoarseFaces_, totalNCoarseFaces_, 0.0)
        );	
		
        skyViewCoeffMatrix_.reset
        (
            new scalarRectangularMatrix(skyViewCoeffSize, totalNCoarseFaces_, 0.0)
        );	

        sunViewCoeffMatrix_.reset
        (
            new scalarRectangularMatrix(sunViewCoeffSize, totalNCoarseFaces_, 0.0)
        );

        Info<< "Insert elements in the matrix..." << endl;

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            insertMatrixElements
            (
                globalNumbering,
                procI,
                globalFaceFacesProc[procI],
                F[procI],
                Fmatrix_()
            );
        }
		
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            insertRectangularMatrixElements
            (
                globalNumbering,
                procI,
                sunskyMap,
                globalFaceFacesProc[procI],
                skyViewCoeff[procI],
                skyViewCoeffMatrix_()
            );
        }

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            insertRectangularMatrixElements
            (
                globalNumbering,
                procI,
                sunskyMap,
                globalFaceFacesProc[procI],
                sunViewCoeff[procI],
                sunViewCoeffMatrix_()
            );
        }			

        bool smoothing = readBool(coeffs_.lookup("smoothing"));
        if (smoothing)
        {
            Info<< "Smoothing the matrix..." << endl;

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                scalar sumF = 0.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    sumF += Fmatrix_()[i][j];
                }
                scalar delta = sumF - 1.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    Fmatrix_()[i][j] *= (1.0 - delta/sumF);
                }
            }
        }

        constEmissivity_ = readBool(coeffs_.lookup("constantEmissivity"));
        if (constEmissivity_)
        {
            CLU_.reset
            (
                new scalarSquareMatrix
                (
                    totalNCoarseFaces_,
                    totalNCoarseFaces_,
                    0.0
                )
            );

            pivotIndices_.setSize(CLU_().n());
        }
	
        timestepsInADay_ = readLabel(coeffs_.lookup("timestepsInADay"));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::directAndDiffuse::directAndDiffuse(const volScalarField& T)
:
    solarLoadModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    Qs_
    (
        IOobject
        (
            "Qs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
	skyViewCoeffMatrix_(),
	sunViewCoeffMatrix_(),	
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),	
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),	
    constEmissivity_(false),
    timestepsInADay_(24),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


Foam::solarLoad::directAndDiffuse::directAndDiffuse
(
    const dictionary& dict,
    const volScalarField& T
)
:
    solarLoadModel(typeName, dict, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    Qs_
    (
        IOobject
        (
            "Qs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
	skyViewCoeffMatrix_(),
	sunViewCoeffMatrix_(),	
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),	
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),	
    constEmissivity_(false),
    timestepsInADay_(24),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarLoad::directAndDiffuse::~directAndDiffuse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solarLoad::directAndDiffuse::read()
{
    if (solarLoadModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::solarLoad::directAndDiffuse::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(procI, faceI);
        forAll(globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }

    //Info << "Fmatrix: " << Fmatrix << endl;
}

void Foam::solarLoad::directAndDiffuse::insertRectangularMatrixElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& sunskyMap,
    const labelListList& globalFaceFaces,
    const scalarListList& skysunViewCoeffs,
    scalarRectangularMatrix& skysunViewCoeffMatrix
)
{
    forAll(skysunViewCoeffs, vectorId)
    {
        const scalarList& vf = skysunViewCoeffs[vectorId];

        //const labelList& globalFaces = globalFaceFaces[vectorId];
        //label globalI = globalNumbering.toGlobal(procI, vectorId);

        //Info << "globalFaces: " << globalFaces << endl; 

        forAll(vf, faceI)
        {        
            skysunViewCoeffMatrix[vectorId][sunskyMap[procI][faceI]] = vf[faceI];
        }
    }

   //Info << "skysunViewCoeffMatrix: " << skysunViewCoeffMatrix << endl;
}

void Foam::solarLoad::directAndDiffuse::calculate()
{
    // Store previous iteration
    Qs_.storePrevIter();

    scalarField compactCoarseT(map_->constructSize(), 0.0);
    scalarField compactCoarseE(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Fill local averaged(T), emissivity(E) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseTave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseEave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    forAll(selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& Tp = T_.boundaryField()[patchID];
        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        fvPatchScalarField& QsPatch = Qs_.boundaryField()[patchID];

        solarLoadViewFactorFixedValueFvPatchScalarField& Qsp =
            refCast
            <
                solarLoadViewFactorFixedValueFvPatchScalarField
            >(QsPatch);

        const scalarList eb = Qsp.albedo();

        const scalarList& Hoi = Qsp.Qso();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID]; 
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID]; 

        scalarList Tave(pp.size(), 0.0);
        scalarList Eave(Tave.size(), 0.0);
        scalarList Hoiave(Tave.size(), 0.0);

        if (pp.size() > 0)
        {
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom) + 1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            forAll(coarseToFine, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                UIndirectList<scalar> fineSf
                (
                    sf,
                    fineFaces
                );
                scalar area = sum(fineSf());
                // Temperature, emissivity and external flux area weighting
                forAll(fineFaces, j)
                {
                    label faceI = fineFaces[j];
                    Tave[coarseI] += (Tp[faceI]*sf[faceI])/area;
                    Eave[coarseI] += (eb[faceI]*sf[faceI])/area;
                    Hoiave[coarseI] += (Hoi[faceI]*sf[faceI])/area;
                }
            }
        }

        localCoarseTave.append(Tave);
        localCoarseEave.append(Eave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseT,nLocalCoarseFaces_).assign(localCoarseTave);
    SubList<scalar>(compactCoarseE,nLocalCoarseFaces_).assign(localCoarseEave);
    SubList<scalar>
        (compactCoarseHo,nLocalCoarseFaces_).assign(localCoarseHoave);

    // Distribute data
    map_->distribute(compactCoarseT);
    map_->distribute(compactCoarseE);
    map_->distribute(compactCoarseHo);

    // Distribute local global ID
    labelList compactGlobalIds(map_->constructSize(), 0.0);

    labelList localGlobalIds(nLocalCoarseFaces_);

    for(label k = 0; k < nLocalCoarseFaces_; k++)
    {
        localGlobalIds[k] = globalNumbering.toGlobal(Pstream::myProcNo(), k);
    }

    SubList<label>
    (
        compactGlobalIds,
        nLocalCoarseFaces_
    ).assign(localGlobalIds);

    map_->distribute(compactGlobalIds);

    // Create global size vectors
    scalarField T(totalNCoarseFaces_, 0.0);
    scalarField E(totalNCoarseFaces_, 0.0);
    scalarField QsExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll(compactCoarseT, i)
    {
        T[compactGlobalIds[i]] = compactCoarseT[i];
        E[compactGlobalIds[i]] = compactCoarseE[i];
        QsExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(T, maxEqOp<scalar>());
    Pstream::listCombineGather(E, maxEqOp<scalar>());
    Pstream::listCombineGather(QsExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(T);
    Pstream::listCombineScatter(E);
    Pstream::listCombineScatter(QsExt);

    // Net solarLoad
    scalarField q(totalNCoarseFaces_, 0.0);

    if (Pstream::master())
    {
        // Variable emissivity
        if (!constEmissivity_)
        {
            scalarSquareMatrix C(totalNCoarseFaces_, totalNCoarseFaces_, 0.0);

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    scalar invEj = 1.0/E[j];
                    scalar sigmaT4 =
                        physicoChemical::sigma.value()*pow(T[j], 4.0);

                    if (i==j)
                    {
                        C[i][j] = invEj - (invEj - 1.0)*Fmatrix_()[i][j];
                        q[i] += (Fmatrix_()[i][j] - 1.0)*sigmaT4 - QsExt[j];
                    }
                    else
                    {
                        C[i][j] = (1.0 - invEj)*Fmatrix_()[i][j];
                        q[i] += Fmatrix_()[i][j]*sigmaT4 - QsExt[j];
                    }

                }
            }

            Info<< "\nSolving view factor equations..." << endl;
            // Negative coming into the fluid
            LUsolve(C, q);
        }
        else //Constant emissivity
        {
            //dimensionedScalar S(coeffs_.lookup("S"));           
            //dimensionedScalar D(coeffs_.lookup("D"));   
            
            //Info << "sunViewCoeffList_: " << sunViewCoeffList_() << endl;
            //Info << "skyViewCoeffList_: " << skyViewCoeffList_() << endl;
            Time& time = const_cast<Time&>(mesh_.time());
            Info << "time.value(): " << time.value() << endl; Info << "mesh_.time(): " << mesh_.time().value() << endl;
            label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 ); Info << "1timestep: " << timestep; timestep = timestep%timestepsInADay_; Info << ", 2timestep: " << timestep << endl;
            //Info << "timestep: " << timestep << endl;
            //Info << "sunViewCoeffList_()[timestep][3]: " << sunViewCoeffList_()[timestep][3] << endl;

            // Initial iter calculates CLU and chaches it
            if (iterCounter_ == 0)
            {
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    { 
                        //scalar invEj = 1/E[j];
                        if (i==j)
                        {
							//CLU_()[i][j] = (1/(1-E[j]));//+(E[j]/(1-E[j]))*Fmatrix_()[i][j];
                            CLU_()[i][j] = (1/(1-E[j]))-(E[j]/(1-E[j]))*Fmatrix_()[i][j];
                        }
                        else
                        {
                            //CLU_()[i][j] = (E[j]/(1-E[j]))*Fmatrix_()[i][j];
							CLU_()[i][j] = -(E[j]/(1-E[j]))*Fmatrix_()[i][j];
                        }
                    }
                }
                Info<< "\nDecomposing C matrix..." << endl;
                LUDecompose(CLU_(), pivotIndices_);
            }
			
            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
					scalar Id = (skyViewCoeffMatrix_()[timestep][j] + sunViewCoeffMatrix_()[timestep][j]);
                    if (i==j)
                    {
                        //q[i] += (Fmatrix_()[i][j] + 1.0)*Id - QsExt[j];
                        q[i] += (Fmatrix_()[i][j] - 1.0)*(-Id) - QsExt[j];
                    }
                    else
                    {
                        //q[i] += (Fmatrix_()[i][j])*Id - QsExt[j];
                        q[i] += (Fmatrix_()[i][j])*(0) - QsExt[j];
                    }
                }
            }

            Info<< "\nLU Back substitute C matrix.." << endl;
            LUBacksubstitute(CLU_(), pivotIndices_, q);
            iterCounter_ ++;
        }
    }

    // Scatter q and fill Qs
    Pstream::listCombineScatter(q);
    Pstream::listCombineGather(q, maxEqOp<scalar>());


    label globCoarseId = 0;
    forAll(selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& Qsp = Qs_.boundaryField()[patchID];
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_.patchFaceMap()[patchID];

            scalar heatFlux = 0.0;
            forAll(coarseToFine, coarseI)
            {
                label globalCoarse =
                    globalNumbering.toGlobal(Pstream::myProcNo(), globCoarseId);
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                forAll(fineFaces, k)
                {
                    label faceI = fineFaces[k];

                    Qsp[faceI] = q[globalCoarse];
                    heatFlux += Qsp[faceI]*sf[faceI];
                }
                globCoarseId ++;
            }
        }
    }

    if (debug)
    {
        forAll(Qs_.boundaryField(), patchID)
        {
            const scalarField& Qsp = Qs_.boundaryField()[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            scalar heatFlux = gSum(Qsp*magSf);
            Info<< "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax Qs if necessary
    Qs_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::solarLoad::directAndDiffuse::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::solarLoad::directAndDiffuse::Ru() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}

// ************************************************************************* //
