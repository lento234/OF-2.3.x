
    labelList viewFactorsPatches(patches.size());
    label count = 0;

    forAll(LADb, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& LADpI = LADb[patchI];

    		if ((isA<wallFvPatch>(mesh.boundary()[patchI])) && (pp.size() > 0))
        {
            viewFactorsPatches[count] = LADpI.patch().index();
            count ++;
        }
    }
    viewFactorsPatches.resize(count--);

    //

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

            //point sf = coarseSf[faceI];
            localCoarseCf.append(cf);
            //localCoarseSf.append(sf);
            //localAgg.append(agglomI);
        }
    }

    //Info << "localCoarseCf " <<  localCoarseCf << endl;
