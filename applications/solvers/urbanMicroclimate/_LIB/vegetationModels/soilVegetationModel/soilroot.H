namespace Foam
{

//- Calculate soil-root system conductance (s/m)
void soilVegetationModel::calc_gsr(const volScalarField& Kl)
{
    // Calculate conductivies
    scalar K, ks, kr, meshDL;
    forAll(RAD_, cellI)
    {
        if (RAD_[cellI] > minThreshold)
        {
            //- Hydraulic conductivity [m/s]
            K = Kl[cellI] * gabs_.value();

            // Soil-root interface conductance [1/s]
            ks = alpha_.value() * K * RAD_[cellI];

            // Characterstic cell size
            meshDL = pow(meshSoil_.V()[cellI],1.0/3.0);

            // Root system condutance (1/s)
            kr = RAD_[cellI] * meshDL / beta_.value();

            // Effective soil-root conductivity (s/m)
            gsr_[cellI] = (ks * kr) / ((ks + kr) * gabs_.value());
        }
    }
    gsr_.correctBoundaryConditions();

    //Info << "Debug ::: min = " << gMin(gsr_) << ", max = " << gMax(gsr_) << endl;
}

//- Calculate root water potential (Pa)
scalar soilVegetationModel::calc_psi_R(const scalar& psi_L, const scalar& gx)
{
    return psi_L + E_.value() / gx;
}

//- Objective function for SPAC minimization problem
scalar soilVegetationModel::SPAC_objective_fx(const scalar& psi_L, const volScalarField& psi_S)
{
    
    //- Xylem condutance [s/m]
    gx_.value() = calc_gx(psi_L);

    //- Calculate root water potential (Pa)
    psi_R_.value() = calc_psi_R(psi_L, gx_.value());

    //- Calculate root update (kg/m2s)
    gv_root_ = - gsr_ * (psi_S - psi_R_);
    
    gv_root_.correctBoundaryConditions();

    //- Calculate net root uptake (kg/s)
    const dimensionedScalar Qr = fvc::domainIntegrate( RAD_ * gv_root_);

    //- Error (kg/s)
    scalar error = E_.value() + Qr.value();

    // Info << "E = " << E_.value()
    //      << ", Qr =  " << Qr.value()
    //      << ", error = " << error << endl;

    return error;
}


// solve Soil-Plant-Atmosphere continuum (WIP)
void soilVegetationModel::solve_SPAC(const volScalarField& Kl, const volScalarField& psi_S)
{

    //- Calculate effective soil-root system conductance (s/m) - gsr_
    calc_gsr(Kl);

    // Info << "gsr: min  = " << gMin(gsr_)
    //      << ", max = " << gMax(gsr_)
    //      << endl;
         
    scalar error, relError, xn, fxnm1, fxnm2;
    scalar xnm1 = 0;
    scalar xnm2 = -10e6;
    
    scalar reltol = 1e-12; //#TODO#: should be in solutiondict
    label maxiter = 200;   //#TODO#: should be in solutiondict

    // Solve the objective function
    fxnm2 = SPAC_objective_fx(xnm2, psi_S); // obtain the mismatch

    // Solve for leaf water potential
    label i;
    for (i=0; i < maxiter; i++)
    {
        fxnm1 = SPAC_objective_fx(xnm1, psi_S);

        // Secant method
        if (mag(fxnm1-fxnm2) < SMALL)
        { 
            Info << "[WARNING] :: oh oh.... stranger danger" << endl;
            fxnm1 = SPAC_objective_fx(xnm1*0.5, psi_S);
        } 
        xn =  xnm1 - fxnm1 * (xnm1 - xnm2) / (fxnm1 - fxnm2);

        // error = mag(fxnm1);
        // relError = mag(fxnm1) / mag( E_.value() );
        
        // relative error
        error = mag(xn-xnm1);
        relError = error / mag(xn);
        
        // Info << "secant iter i = " << i 
        //      << ", psi_L = " << xn/1e6 
        //      << " MPa, error = " << error
        //      << ", rel. error = " << relError << endl;

        // Check convergence
        if ( relError < reltol)
        {
            // Info << "converged ! error = " << relError << endl;
            break;
        }
        else
        {
            xnm2 = xnm1;
            xnm1 = xn;
            fxnm2 = fxnm1;
        }
        
    }

    // Update psi_L
    psi_L_.value() = xnm1;

    
    //#TODO#: Remove, if okay
    if ( (psi_L_.value() > 0) || (psi_L_.value() < -10e6) )
    {
        Info << "                       [WARNING] "
             << " psi_L outside bound, "
             << "psi_L = " << psi_L_.value() 
             << endl;
    }
    if ( (psi_R_.value() > 0) )
    {
        Info << "                       [WARNING] "
             << " psi_R positive, "
             << "psi_R = " << psi_R_.value() 
             << endl;
    }
    if ( (psi_R_.value() < psi_L_.value()) )
    {
        Info << "                       [WARNING] "
             << " psi_R < psi_L (reverse), "
             << "psi_L = " << psi_L_.value() 
             << "psi_R = " << psi_R_.value() 
             << endl;
    }



    Info << "VEGETATION: [SPAC] : psi_L = " << round(1e4 * psi_L_.value() / 1e6) / 1e4
         << " MPa, psi_R = " << round( 1e4 * psi_R_.value() / 1e6 ) / 1e4
         << " MPa, psi_S = " << round( 1e4 * average(psi_S).value() / 1e6 ) / 1e4
         << " MPa, Final rel. residual = " << relError
         << ", No Iterations " << i << endl;
    
}

}