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
}

//- Calculate soil xylem condutance (1/s)
scalar soilVegetationModel::calc_gx(const scalar& psi_L)
{
    // Xylem condutance (1/s)
    scalar gx = gx_max_.value() * exp( - pow(-psi_L/d_.value(),cx_) );
    
    // Effective xylem condutance (sm)
    return Ax_.value() * gx / gabs_.value();
}

//- Calculate marginal water use efficiency
void soilVegetationModel::calc_marginalWUE()
{
    // Average leaf water potential of the last 24 hours
    scalar psi_L24_avg = average(psi_L24_);

    Info << "Vegetation : [MWUE]        :: psi_L24_avg =  " << psi_L24_avg  << endl;

    // Calculate Marginal Water use efficiency (lamba) 
    lambda_.value() = lambda_max_.value() * (ca_.value() / ca_star_.value()) * exp ( - betaL_.value() * pow(psi_L24_avg - psi_Lmax_.value(), 2));
}

/*
//- Calculate soil xylem condutance (1/s) (WIP)
scalar soilVegetationModel::calc_psi_R(const scalar& psi_L, const scalar& gx)
{
    return psi_L + E_.value() / gx;
}

//- Solve SPAC error (WIP)
scalar soilVegetationModel::SPAC_fx(const scalar& psi_L)
{
    
    //- Xylem condutance [s/m]
    gx_.value() = calc_gx(psi_L);

    //- Calculate root water potential (Pa)
    psi_R_.value() = calc_psi_R(psi_L, gx_.value());

    //- Calculate net root update
    Qr_ = fvc::domainIntegrate( -gsr_ * (psi_S - psi_R_) * RAD_ );

    //- Error
    scalar error = E_.value() + Qr_.value();

    return error;
}


// solve Soil-Plant-Atmosphere continuum (WIP)
void soilVegetationModel::solve_SPAC(const volScalarField& pc, const volScalarField& Kl)
{
    //- Net transpiration rate (kg/s) 
    E_.value() = fvc::domainIntegrate(LAD_ * gv_leaf_);


    //- Calculate effective soil-root system conductance (s/m) - gsr_
    calc_gsr(Kl);

    //- Calculate soil water potential
    volScalarField& psi_S = pc ; // + g & meshSoil_.C();

    scalar psi_L;
    scalar error;

    scanar xn;
    scalar xnm1 = 0;
    scalar xnm2 = -10e6;
    scalar fxnm1, fxnm2;

    scalar reltol = 1e-12;
    label maxiter = 100;

    fxnm2 = SPAC_fx(xnm2);

    // Solve for leaf water potential
    for (label i=0; i < maxiter; i++)
    {
        fxnm1 = SPAC_fx(xnm1);

        // Secant
        xn =  xnm1 - fxnm1 * (xnm1 - xnm2) / (fxnm1 - fxnm2);

        // Check convergence
        if ( (abs(xn-xnm1) / abs(xn)) < reltol)
        {
            break;
        }
        else
        {
            xnm2 = xnm1;
            xnm1 = xn;
            fxnm2 = fxnm1;
        }
    }

    psi_L_.value() = xnm1;

    // Update psi_L24
    for (label i=psi_L24_.size()-1; i>0; i--)
    {  
        psi_L24_[i-1] = psi_L24_[i];
    }
    psi_L24_[psi_L24_.size()-1] = psi_L_.value();

}

*/

}