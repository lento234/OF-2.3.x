// All rays expressed as start face (local) index end end face (global)
// Pre-size by assuming a certain percentage is visible.

// Maximum length for dynamicList
const label maxDynListLength = 10000;

//for (label procI = 0; procI < Pstream::nProcs(); procI++)
//{
    // Shoot rays from me to procI. Note that even if processor has
    // 0 faces we still need to call findLine to keep calls synced.

    DynamicField<point> start(coarseMesh.nFaces());
    DynamicField<vector> startArea(coarseMesh.nFaces());
    DynamicField<point> end(start.size());
    DynamicList<label> startIndex(start.size());
    DynamicList<label> endIndex(start.size());

    const pointField& myFc = localCoarseCf_[Pstream::myProcNo()];
    const vectorField& myArea = localCoarseSf_[Pstream::myProcNo()];

    //const pointField& remoteArea = remoteCoarseSf[procI];
    const pointField& remoteFc = remoteCoarseCf_[Pstream::myProcNo()];

    label i = 0;
    do
    {
        for (; i < myFc.size(); i++)
        {
            const point& fc = myFc[i];
            const vector& area = myArea[i];
            const point& remFc = remoteFc[i];
            
            const vector& d = remFc - fc;
           
            start.append(fc + 0.0001*d);
            startIndex.append(i);
            startArea.append(area);
            end.append(fc + 0.9999*d);
            label globalI = globalNumbering.toGlobal(Pstream::myProcNo(), i);
            endIndex.append(globalI);
            if (startIndex.size() > maxDynListLength)
            {
                break;
            }
        }
        
        List<pointIndexHit> hitInfo(startIndex.size());
        surfacesMesh.findLine(start, end, hitInfo);
        
        Info << "start: " << start << endl;
        //Info << "startArea: " << startArea << endl;
        Info << "end: " << end << endl;
        
        forAll (hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit()) // if the ray did not hit any surface between start and end positions...
            {
                scalar cosPhiTest = (startArea[rayI] & sunPos)/(mag(startArea[rayI])*mag(sunPos) + SMALL);
                if (cosPhiTest < 0) // ... and if the surface is looking towards the sun
                {
                    rayStartFace.append(startIndex[rayI]);
                    rayEndFace.append(endIndex[rayI]);
                }
            }
        }

        start.clear();
        startIndex.clear();
        end.clear();
        endIndex.clear();

    }while (returnReduce(i < myFc.size(), orOp<bool>()));

//}
