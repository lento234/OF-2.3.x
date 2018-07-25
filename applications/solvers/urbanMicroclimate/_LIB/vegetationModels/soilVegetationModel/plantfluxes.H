namespace Foam
{

//- Calculate soil xylem condutance (1/s)
scalar soilVegetationModel::calc_gx(const scalar& psi_L)
{
    // Xylem condutance (1/s)
    scalar gx = gx_max_.value() * exp( - pow(-psi_L/d_.value(),cx_) );
    
    // Effective xylem condutance (sm)
    return Ax_.value() * gx / gabs_.value();
}


void soilVegetationModel::update_psiL24()
{
    for (label i=0; i<psi_L24_.size(); i++)
    {  
        psi_L24_[i] = psi_L24_[i+1];
    }
    psi_L24_[psi_L24_.size()-1] = psi_L_.value();
}

//- Calculate marginal water use efficiency
void soilVegetationModel::calc_marginalWUE()
{
    // Update psi_L24
    update_psiL24();
    
    // Average leaf water potential of the last 24 hours
    scalar psi_L24_avg = average(psi_L24_);

    Info << "Vegetation : [MWUE]        :: psi_L24_avg =  " << psi_L24_avg  << endl;

    // Calculate Marginal Water use efficiency (lamba) 
    lambda_.value() = lambda_max_.value() * (ca_.value() / ca_star_.value()) * exp ( - betaL_.value() * pow(psi_L24_avg - psi_Lmax_.value(), 2));
    //lambda_.value() = lambda_max_.value() * (ca_.value() / ca_star_.value()) * exp ( - betaL_.value() * psi_L24_avg ); * betaL (1/MPa)
}

} // FOAM
