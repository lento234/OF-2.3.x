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
            end.append(fc + 1.0001*d);
            label globalI = globalNumbering.toGlobal(Pstream::myProcNo(), i);
            endIndex.append(globalI);
            if (startIndex.size() > maxDynListLength)
            {
                break;
            }
        }
        
        List<pointIndexHit> hitInfo(startIndex.size());
        List<pointIndexHit> hitInfoCyclic(startIndex.size());

	    bool cyclicTest = true;   
        do
        {	 
		cyclicTest = false;
      
		surfacesMesh.findLine(start, end, hitInfo);
		cyclicMesh.findLine(start, end, hitInfoCyclic); 
		
        vectorField normals(startIndex.size());
		cyclicMesh.getNormal(hitInfoCyclic,normals);
		
        Info << "start: " << start << endl;

        forAll (start, rayI)
		{
		    Info << "rayI: " << rayI << start[rayI] << end[rayI] << endl;
			
		    if (!hitInfo[rayI].hit()) // if the ray did not hit any surface between start and end positions...
		    {
		        scalar cosPhiTest = (startArea[rayI] & sunPos)/(mag(startArea[rayI])*mag(sunPos) + SMALL);
		        if (cosPhiTest < 0) // ... and if the surface is looking towards the sun
		        {
		            if (hitInfoCyclic[rayI].hit()) // if the ray hit cyclic boundaries
		            {
		                //Info << "Startpoint was here: " << start[rayI] << endl;
		                //Info << "Endpoint was here: " << end[rayI] << endl;

                        vector normal = normals[rayI];
       	                vector diff = max_ - min_;
	                    vector diffDotNormal = ((diff & normal) * cmptMag(normal));

                        point newFc = hitInfoCyclic[rayI].hitPoint() - diffDotNormal;

		                scalar i1 = 0; scalar i2 = 0; scalar i3 = 0;
		                if (sunPos.x() > 0.0)
		                {
		                    i1 = (max_.x() - newFc.x())/sunPos.x();
		                } 
		                else if (sunPos.x() < 0.0)
		                {
		                    i1 = (min_.x() - newFc.x())/sunPos.x();
		                } 
		                else {i1 = VGREAT;}

		                if (sunPos.y() > 0.0)
		                {
		                    i2 = (max_.y() - newFc.y())/sunPos.y();
		                } 
		                else if (sunPos.y() < 0.0)
		                {
		                    i2 = (min_.y() - newFc.y())/sunPos.y();
		                }
		                else{i2 = VGREAT;}

		                if (sunPos.z() > 0.0)
		                {
		                    i3 = (max_.z() - newFc.z())/sunPos.z();
		                } 
		                else if (sunPos.z() < 0.0)
		                {
		                    i3 = (min_.z() - newFc.z())/sunPos.z();
		                }
		                else{i3 = VGREAT;}
                            //Info << "newFc: " << newFc << endl;
		                scalar k = min(i1, min(i2, i3));
                        point newRemFc = k*point(sunPos.x(),sunPos.y(),sunPos.z())+newFc;
                        vector newD = newRemFc - newFc; 
		                point newStart = newFc + 0.0001*newD;
                        point newEnd = newFc + 1.0001*newD;
                            Info << "New Startpoint is here: " << newStart << endl;
		                    //Info << "New Endpoint is here: " << newEnd << endl;
		                start[rayI] = newStart;
		                end[rayI] = newEnd;

		                cyclicTest = true;
		            }
		            else
		            {
		                rayStartFace.append(startIndex[rayI]);
		                rayEndFace.append(endIndex[rayI]);
                        Info << "Clear" << endl;
		            }
		        }
		        else
		        {
		            Info << "Not looking towards the sun" << endl;
		        }
            }
		    else
		    {
		        Info << "Hit a wall" << endl;
		    }
					
		}

		if(cyclicTest == true)
		{
		    rayStartFace.clear();
		    rayEndFace.clear();                 
		}
	}while (returnReduce(cyclicTest, orOp<bool>()));               

        start.clear();
        startIndex.clear();
        end.clear();
        endIndex.clear();
        hitInfo.clear();
        hitInfoCyclic.clear();
    }while (returnReduce(i < myFc.size(), orOp<bool>()));

//}
