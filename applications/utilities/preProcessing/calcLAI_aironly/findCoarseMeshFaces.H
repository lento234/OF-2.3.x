
    labelList viewFactorsPatches(patches.size());
    labelList howManyCoarseFacesPerPatch(patches.size());
    label count = 0;
    label countAll = 0;

    //forAll(LADb, patchI)
    forAll(Qrb, patchI)
    {
        const polyPatch& pp = patches[patchI];
        //const fvPatchScalarField& LADpI = LADb[patchI];
        const fvPatchScalarField& QrpI = Qrb[patchI];

    		if ((isA<wallFvPatch>(mesh.boundary()[patchI])) && (pp.size() > 0))
        {
            //viewFactorsPatches[count] = LADpI.patch().index();
            viewFactorsPatches[count] = QrpI.patch().index();
            count ++;
            howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchI].size();
        }
        else if ((isA<fixedValueFvPatchScalarField>(QrpI)) && (pp.size() > 0))
		    {
            howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchI].size();
        }
        else
        {
            howManyCoarseFacesPerPatch[countAll] = 0;
        }
        countAll ++;
    }

    viewFactorsPatches.resize(count--);

    DynamicList<point> localCoarseCf(nCoarseFaces);

    forAll (viewFactorsPatches, i)
    {
        const label patchID = viewFactorsPatches[i];

        const polyPatch& pp = patches[patchID];
        const labelList& agglom = finalAgglom[patchID];

        label nAgglom = max(agglom)+1;
        labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
        const labelList& coarsePatchFace = coarseMesh.patchFaceMap()[patchID];

        const pointField& coarseCf = coarseMesh.Cf().boundaryField()[patchID];
        //const pointField& coarseSf = coarseMesh.Sf().boundaryField()[patchID];

        labelHashSet includePatches;
        includePatches.insert(patchID);

        forAll(coarseCf, faceI)
        {
            point cf = coarseCf[faceI];

            const label coarseFaceI = coarsePatchFace[faceI];
            const labelList& fineFaces = coarseToFine[coarseFaceI];
            // const label agglomI =
            //     agglom[fineFaces[0]] + coarsePatches[patchID].start();

            // Construct single face
            uindirectPrimitivePatch upp
            (
                UIndirectList<face>(pp, fineFaces),
                pp.points()
            );

            List<point> availablePoints
            (
                upp.faceCentres().size()
              + upp.localPoints().size()
            );

            SubList<point>
            (
                availablePoints,
                upp.faceCentres().size()
            ).assign(upp.faceCentres());

            SubList<point>
            (
                availablePoints,
                upp.localPoints().size(),
                upp.faceCentres().size()
            ).assign(upp.localPoints());

            point cfo = cf;
            scalar dist = GREAT;
            forAll(availablePoints, iPoint)
            {
                point cfFine = availablePoints[iPoint];
                if(mag(cfFine-cfo) < dist)
                {
                    dist = mag(cfFine-cfo);
                    cf = cfFine;
                }
            }
            localCoarseCf.append(cf);
        }
    }

    //Info << "localCoarseCf " <<  localCoarseCf << endl;
