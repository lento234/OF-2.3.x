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

}