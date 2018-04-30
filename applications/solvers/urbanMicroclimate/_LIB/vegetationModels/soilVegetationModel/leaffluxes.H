
namespace Foam
{

void soilVegetationModel::solve_assimilation(const scalar& Tl, const scalar& VPD, const scalar& c, const label& cellI)
{
    // Aerodynamic resistance
    scalar ga = ga_[cellI];

    // Calculating rubisco limited coefficients (A_C)

        // Coefficients
        scalar Vcmax = Vcmax25_.value() * exp(0.088*(Tl-273.15-25.0))/(1.0 + exp(0.29*(Tl-273.15-41.0)));
        scalar Kc = Kc25_.value() * exp(gammac_*(Tl-273.15-25.0));
        scalar Ko = Ko25_.value() * exp(gammao_*(Tl-273.15-25.0));

        // Dark respiration
        scalar Rd25 = 0.015 * Vcmax25_.value();
        scalar Rd = Rd25 * exp(0.069*(Tl - 273.15 - 25.0)) / (1.0 + exp(1.3 * (Tl-273.15 - 55.0)));

        scalar a1 = Vcmax;
        scalar a2 = Kc * (1.0 + cao_.value() / Ko);
        scalar ccp = (Kc/(2*Ko)) * cao_.value() * (ko_.value()/kco_.value());

        // Condutance due to limited Ac
        scalar gs_C = calc_gs(a1, a2, VPD, c);
        scalar gc_eff_C = calc_gc_eff(gs_C, ga);

        // Intercelluar CO2 due to limited Ac
        scalar ci_C = calc_ci(a1, a2, ccp, gc_eff_C, c, Rd);

        // Assimilaion rate Ac
        scalar A_C = calc_An(gc_eff_C, ci_C, c);

    // Calculate light limiting coefficients (A_E)

        // Quantum flux density
        scalar Qp = rPAR_ * qrsw_leaf_[cellI] / E_PAR_.value();
        

        a1 = gammaPAR_*Qp;
        a2 = 2.0*ccp;

        // Conductance due to limited Ae
        scalar gs_E = calc_gs(a1, a2, VPD, c);
        scalar gc_eff_E = calc_gc_eff(gs_E, ga);

        // Intercelluar CO2 due to limited Ac
        scalar ci_E = calc_ci(a1, a2, ccp, gc_eff_E, c, Rd);

        // Assimilate rate AE
        scalar A_E = calc_An(gc_eff_E, ci_E, c);
        
    // Minimum assimilation
    if (mag(A_E) < mag(A_C))
    {
        // Light limited 
        gs_[cellI] = gs_E;
        gc_eff_[cellI] = gc_eff_E;
        ci_[cellI] = ci_E;
        An_[cellI] = A_E;
    }
    else
    {
        // Rubisco limited 
        gs_[cellI] = gs_C;
        gc_eff_[cellI] = gc_eff_C;
        ci_[cellI] = ci_C;
        An_[cellI] = A_C;       
    }

    // Calculate effective vapour conductance (mol / m2 s)
    gv_eff_[cellI] = calc_gv_eff(gc_eff_[cellI]);

}


//- Calculate Assimilation rate (mol / m2 s)
scalar soilVegetationModel::calc_An(const scalar& gc_eff, const scalar& ci, const scalar& c)
{
    return gc_eff * (ci - c);
}

//- Calculate the intercellular CO2 concentration
scalar soilVegetationModel::calc_ci(const scalar& a1, const scalar& a2, const scalar& ccp, const scalar& gc_eff, const scalar& c, const scalar& Rd)
{
    //scalar qb = (a1 / gc_eff) + a2 - c;
    //scalar qc = - (a1 * ccp / gc_eff) - c*a2;
    scalar qb = (a1 / gc_eff) + a2 - c - Rd/gc_eff; // added dark respiration
    scalar qc = - (a1 * ccp / gc_eff) - c*a2 - a2*Rd/gc_eff; // added dark respiration (mitochondrial)

    // Calculate intercellular CO2 concentration
    scalar ci = solve_quadratic(qb, qc);

    return ci;
}

//- Calculate aerodynamic condutance (mol/m2 s)
void soilVegetationModel::calc_ga(const volVectorField& U, const volScalarField& rho)
{
    // Magnitude of velocity
    //volScalarField magU("magU", mag(U));
    scalar magU;
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > minThreshold)
        {
            magU = mag(U[cellI]);
            ga_[cellI] = (rho[cellI]/(Mco2_.value()*C_.value())) * pow(magU / l_.value(), 0.5); 
        }
    }

    ga_.correctBoundaryConditions();

}

//- Calculate effective CO2 condutance (mol / m2 s)
scalar soilVegetationModel::calc_gc_eff(const scalar& gs, const scalar& ga)
{
    return (gs * ga) / (gs + ga);
}

//- Calculate stomatal condutance (mol / m2 s)
scalar soilVegetationModel::calc_gs(const scalar& a1, const scalar& a2, const scalar& VPD, const scalar& c)
{
    return a1/(a2 + s_*c + SMALL) * ( -1.0 + sqrt(c/(a_*lambda_.value()*VPD + SMALL)) ) + gsn_.value();
}

//- Calculate effective vapour conductance (mol / m2 s)
scalar soilVegetationModel::calc_gv_eff(const scalar& gc_eff)
{
    return a_ * gc_eff;
}

//- Solve the quadratic equation
scalar soilVegetationModel::solve_quadratic(const scalar& qb,const scalar& qc)
{
    scalar b2m4c = pow(qb,2) - 4*qc;
    
    if (b2m4c < 0.0)
        FatalError << "Root is complex: " << b2m4c << abort(FatalError);
    
    scalar x1 = (-qb + sqrt(b2m4c) )/2.0;
    scalar x2 = (-qb - sqrt(b2m4c) )/2.0;

    return max(x1,x2);
}


} // end namepsace Foam

