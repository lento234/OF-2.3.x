{
    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs.last(),timeDirs.size()-1);

    runTime.setDeltaT(storeFluidDeltaT * initialSolidTimestepFactor); //this always resets to the small timestep for solid, could it be better?

    scalar solidInternalTime = 0;//runTime.deltaT().value();
    scalar timeToOutput = storeFluidDeltaT;
    bool timeStepDecrease = false;

    while ( solidInternalTime < storeFluidDeltaT )
    {
        Info << nl << "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaT().value() << endl;
        Info << "solidInternalTime: " << solidInternalTime << endl;
        //pc.storePrevIter();
        Ts.storePrevIter();

        volScalarField Ts_old = Ts;
        //volScalarField ws_old = ws;

        int update_sequence = 1;

        for (int nInternalIter=1; nInternalIter<=nInternalIterMax; nInternalIter++)
        {
            //Moisture transfer////////////
                	// #include "pcEqn.H"
                  //
                  //   //Firstly, test if all pc values are valid
                  //   if (gMax(pc) >= 0)// || gMax(pc.boundaryField()) >= 0)
                  //   {
                  //       Info << "This is going to crash (pc)! Decreasing timestep and reverting fields..." << endl;
                  //       Info << "Error: gMax(pc): " << gMax(pc) << ", gMax(pc.boundaryField()): " << gMax(pc.boundaryField()) << endl;
                  //       //pc.write();
                  //       timeStepDecrease = true;
                  //       pc=pc.prevIter();
                  //       forAll(pc.boundaryField(), patchi)
                  //       {
                  //           forAll(pc.boundaryField()[patchi], facei)
                  //           {
                  //               pc.boundaryField()[patchi][facei] = pc.prevIter().boundaryField()[patchi][facei];
                  //           }
                  //       }
                  //       Ts=Ts.prevIter();
                  //       forAll(Ts.boundaryField(), patchi)
                  //       {
                  //           forAll(Ts.boundaryField()[patchi], facei)
                  //           {
                  //               Ts.boundaryField()[patchi][facei] = Ts.prevIter().boundaryField()[patchi][facei];
                  //           }
                  //       }
                  //       break;
                  //   }
                  //   pc.correctBoundaryConditions();

            ///////////////////////////////

            //Heat transfer////////////////
            	#include "TsEqn.H"

              //Firstly, test if all Ts values are valid
              if (gMin(Ts) <= 0 || gMin(Ts.boundaryField()) <= 0)
              {
                  Info << "This is going to crash (Ts)! Decreasing timestep and reverting fields..." << endl;
                  timeStepDecrease = true;
                  // pc=pc.prevIter();
                  // forAll(pc.boundaryField(), patchi)
                  // {
                  //     forAll(pc.boundaryField()[patchi], facei)
                  //     {
                  //         pc.boundaryField()[patchi][facei] = pc.prevIter().boundaryField()[patchi][facei];
                  //     }
                  // }
                  Ts=Ts.prevIter();
                  forAll(Ts.boundaryField(), patchi)
                  {
                      forAll(Ts.boundaryField()[patchi], facei)
                      {
                          Ts.boundaryField()[patchi][facei] = Ts.prevIter().boundaryField()[patchi][facei];
                      }
                  }
                  break;
              }
              Ts.correctBoundaryConditions();
            ///////////////////////////////


            //Convergence test/////////////

                //update moisture content for convergence test
                update_sequence = 3;
                //pc.dimensions().reset(dimless);
                //ws.dimensions().reset(dimless);
                #include "updatebuildingMaterials.H"
                //ws.dimensions().reset(dimensionSet( 1, -3, 0, 0, 0, 0, 0));
                //pc.dimensions().reset(dimensionSet( 1, -1, -2, 0, 0, 0, 0));

                //convergence test

                //scalar maxChangews = gMax(mag(ws.internalField()-ws_old.internalField()));
                scalar maxChangeTs = gMax(mag(Ts.internalField()-Ts_old.internalField()));
                //if(maxChangews < 0.1 && maxChangeTs < 0.1 && nInternalIter>=2) //force at least 1 internal iteration
                if(maxChangeTs < 0.1 && nInternalIter>=2) //force at least 1 internal iteration
                {
                    //ws_old = ws; ws_old.correctBoundaryConditions();
                    Ts_old = Ts; Ts_old.correctBoundaryConditions();
                    Info << "Total internal iterations: " << nInternalIter << endl;
                    timeStepDecrease = false;
                    break;
                }
                else if (nInternalIter == nInternalIterMax) //nonlinear iteration reached maximum
                {
                    Info
                       << "Nonlinear iteration didn't converge !"
                       << endl
                       //<< "maxChangews: " << maxChangews << endl
                       << "maxChangeT: " << maxChangeTs << endl
                       << "Total internal iterations: " << nInternalIter << endl;
                    timeStepDecrease = true;
                    // pc=pc.prevIter();
                    // forAll(pc.boundaryField(), patchi)
                    // {
                    //     forAll(pc.boundaryField()[patchi], facei)
                    //     {
                    //         pc.boundaryField()[patchi][facei] = pc.prevIter().boundaryField()[patchi][facei];
                    //     }
                    // }
                    //pc.correctBoundaryConditions();
                    Ts=Ts.prevIter();
                    forAll(Ts.boundaryField(), patchi)
                    {
                        forAll(Ts.boundaryField()[patchi], facei)
                        {
                            Ts.boundaryField()[patchi][facei] = Ts.prevIter().boundaryField()[patchi][facei];
                        }
                    }
                    //Ts.correctBoundaryConditions();
                }
                else //not converged nor reached the maximum iteration yet, continue
                {
                    //ws_old = ws; ws_old.correctBoundaryConditions();
                    Ts_old = Ts; Ts_old.correctBoundaryConditions();
                    //Info << "maxChangews: " << maxChangews << ", maxChangeTs: " << maxChangeTs << endl;
                    Info << "maxChangeTs: " << maxChangeTs << endl;
                }
            ///////////////////////////////

        }

        if (timeStepDecrease != true)
        {
            solidInternalTime += runTime.deltaT().value();
            timeToOutput = storeFluidDeltaT - solidInternalTime;
            if (timeToOutput >= 0.0)
            {
                runTime++;
            }
        }

        if (timeToOutput > 0.0)
        {
            //#include "readSolidTimeControls.H"
            #include "setSolidRegionDeltaT.H"

            //update thermal radiation fluxes
            forAll(fluidRegions, i)
            {
                Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
                radiation::radiationModel& rad = radiation[i];
                rad.correct();
            }
        }

        Info << "timeToOutput: " << timeToOutput << endl;
    }

    runTime.setDeltaT(storeFluidDeltaT);
}
