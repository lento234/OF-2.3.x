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

    //- Calculate net root uptake (kg/s)
    const dimensionedScalar Qr_ = fvc::domainIntegrate( RAD_ * gv_root_);

    //- Error (kg/s)
    scalar error = E_.value() + Qr_.value();

    return error;
}


// solve Soil-Plant-Atmosphere continuum (WIP)
void soilVegetationModel::solve_SPAC(const volScalarField& Kl, const volScalarField& psi_S)
{

    //- Calculate effective soil-root system conductance (s/m) - gsr_
    calc_gsr(Kl);

    scalar relError, xn, fxnm1, fxnm2;
    scalar xnm1 = 0;
    scalar xnm2 = -10e6;
    
    scalar reltol = 1e-12; // should be in solutiondict
    label maxiter = 100;   // should be in solutiondict

    // Solve the objective function
    fxnm2 = SPAC_objective_fx(xnm2, psi_S); // obtain the mismatch

    // Solve for leaf water potential
    label i;
    for (i=0; i < maxiter; i++)
    {
        fxnm1 = SPAC_objective_fx(xnm1, psi_S);

        // Secant method
        xn =  xnm1 - fxnm1 * (xnm1 - xnm2) / (fxnm1 - fxnm2);

        // relative error
        relError = mag(xn-xnm1) / mag(xn);
        // Check convergence
        if ( relError < reltol)
        {
            break;
        }
        else
        {
            xnm2 = xnm1;
            xnm1 = xn;
            fxnm2 = fxnm1;
        }
        
        // Info << "secant iter i = " << i 
        //      << ", psi_L = " << xn/1e6 
        //      << " MPa, error = " << relError << endl;
    }

    // Update psi_L
    psi_L_.value() = xnm1;


    Info << "VEGETATION: [SPAC] : psi_L = " << round(1e4 * psi_L_.value() / 1e6) / 1e4
         << " MPa, psi_R = " << round( 1e4 * psi_R_.value() / 1e6 ) / 1e4
         << " MPa, psi_S = " << round( 1e4 * average(psi_S).value() / 1e6 ) / 1e4
         << " MPa, Final rel. residual = " << relError
         << ", No Iterations " << i << endl;
    
}

}