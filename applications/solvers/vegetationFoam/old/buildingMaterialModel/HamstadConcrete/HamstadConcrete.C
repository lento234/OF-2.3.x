/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "HamstadConcrete.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(HamstadConcrete, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        HamstadConcrete,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::HamstadConcrete::HamstadConcrete
(
    const word& name,
    const dictionary& buildingMaterialProperties,
    const word& cellZoneModel
    //volScalarField& h,
    //volScalarField& theta,
    //volScalarField& kr,
    //volScalarField& Ch
)
:
    buildingMaterialModel(name, buildingMaterialProperties, cellZoneModel),// h, theta, kr, Ch),
    HamstadConcreteCoeffs_(buildingMaterialProperties.subDict(typeName + "Coeffs")),
    rho_("rho", dimensionSet(1, -3, 0, 0, 0), HamstadConcreteCoeffs_.lookup("rho")),
    cap_("cap", dimensionSet(0, 2, -2, -1, 0), HamstadConcreteCoeffs_.lookup("cap"))  
    /*Ks_(HamstadConcreteCoeffs_.lookup("Ks")),
    theta_s_(HamstadConcreteCoeffs_.lookup("theta_s")),
    theta_r_(HamstadConcreteCoeffs_.lookup("theta_r")),
    alpha_(HamstadConcreteCoeffs_.lookup("alpha")),
    beta_(HamstadConcreteCoeffs_.lookup("beta")),
    gamma_(HamstadConcreteCoeffs_.lookup("gamma")),
    A_(HamstadConcreteCoeffs_.lookup("A")),
    Ss_(HamstadConcreteCoeffs_.lookup("Ss"))*/
{
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::buildingMaterialModels::HamstadConcrete::read
(
    const dictionary& buildingMaterialProperties
)
{
    buildingMaterialModel::read(buildingMaterialProperties);

    HamstadConcreteCoeffs_ = buildingMaterialProperties.subDict(typeName + "Coeffs");

    HamstadConcreteCoeffs_.lookup("rho") >> rho_.value();
    HamstadConcreteCoeffs_.lookup("cap") >> cap_.value();    

    /*HamstadConcreteCoeffs_.lookup("Ks") >> Ks_;
    HamstadConcreteCoeffs_.lookup("theta_s") >> theta_s_;
    HamstadConcreteCoeffs_.lookup("theta_r") >> theta_r_;
    HamstadConcreteCoeffs_.lookup("alpha") >> alpha_;
    HamstadConcreteCoeffs_.lookup("beta") >> beta_;
    HamstadConcreteCoeffs_.lookup("gamma") >> gamma_;
    HamstadConcreteCoeffs_.lookup("A") >> A_;
    HamstadConcreteCoeffs_.lookup("Ss") >> Ss_;*/

    return true;
}

//- Correct the buildingMaterial moisture content (cell)
void Foam::buildingMaterialModels::HamstadConcrete::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(1); reta[0]=-8.0e-8;
    List<scalar> retn; retn.setSize(1); retn[0]=1.6e0;
    List<scalar> retm; retm.setSize(1); retm[0]=0.375e0;
    List<scalar> retw; retw.setSize(1); retw[0]=1e0;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=0; i++)
    {
        tmp = pow( (reta[i]*pc.internalField()[celli]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.internalField()[celli]); 
    } 
    w.internalField()[celli] = w_tmp*146;   
    Crel.internalField()[celli] = mag( C_tmp*146 );   
}

//- Correct the buildingMaterial moisture content (boundary)
void Foam::buildingMaterialModels::HamstadConcrete::update_w_C_boundary(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label patchi, label patchFacei)
{
    List<scalar> reta; reta.setSize(1); reta[0]=-8.0e-8;
    List<scalar> retn; retn.setSize(1); retn[0]=1.6e0;
    List<scalar> retm; retm.setSize(1); retm[0]=0.375e0;
    List<scalar> retw; retw.setSize(1); retw[0]=1e0;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=0; i++)
    {
        tmp = pow( (reta[i]*pc.boundaryField()[patchi][patchFacei]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.boundaryField()[patchi][patchFacei]); 
    } 
    w.boundaryField()[patchi][patchFacei] = w_tmp*146; 
    Crel.boundaryField()[patchi][patchFacei] = mag( C_tmp*146 );  
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::HamstadConcrete::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar tmp=w.internalField()[celli]-73;
    tmp=-39.2619e0 +0.0704e0*tmp -1.742e-4*pow(tmp,2) -2.7953e-6*pow(tmp,3) -1.1566e-7*pow(tmp,4) +2.5969e-9*pow(tmp,5);
    Krel.internalField()[celli] = pow(10,tmp*0.4342944819e0);

}

//- Correct the buildingMaterial liquid permeability (boundary)
void Foam::buildingMaterialModels::HamstadConcrete::update_Krel_boundary(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label patchi, label patchFacei)
{
    scalar tmp=w.boundaryField()[patchi][patchFacei]-73;
    tmp=-39.2619e0 +0.0704e0*tmp -1.742e-4*pow(tmp,2) -2.7953e-6*pow(tmp,3) -1.1566e-7*pow(tmp,4) +2.5969e-9*pow(tmp,5);
    Krel.boundaryField()[patchi][patchFacei] = pow(10,tmp*0.4342944819e0); 
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::HamstadConcrete::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/146); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*200*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for concrete" [s]
    
    K_v.internalField()[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial vapor permeability (boundary)
void Foam::buildingMaterialModels::HamstadConcrete::update_Kv_boundary(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label patchi, label patchFacei)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.boundaryField()[patchi][patchFacei] - 5.976*Foam::log(T.boundaryField()[patchi][patchFacei])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.boundaryField()[patchi][patchFacei]/(rho_l*R_v*T.boundaryField()[patchi][patchFacei])); // relative humidity [-]
    
    scalar tmp = 1 - (w.boundaryField()[patchi][patchFacei]/146); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.boundaryField()[patchi][patchFacei]*200*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for concrete" [s]
    
    K_v.boundaryField()[patchi][patchFacei] = (delta*p_vsat*relhum)/(rho_l*R_v*T.boundaryField()[patchi][patchFacei]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::HamstadConcrete::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/146); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*200*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for concrete" [s]

    K_pt.internalField()[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (boundary)
void Foam::buildingMaterialModels::HamstadConcrete::update_Kpt_boundary(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label patchi, label patchFacei)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.boundaryField()[patchi][patchFacei] - 5.976*Foam::log(T.boundaryField()[patchi][patchFacei])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.boundaryField()[patchi][patchFacei]*T.boundaryField()[patchi][patchFacei]) - 5.976/T.boundaryField()[patchi][patchFacei]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.boundaryField()[patchi][patchFacei]/(rho_l*R_v*T.boundaryField()[patchi][patchFacei])); // relative humidity [-]
    
    scalar tmp = 1 - (w.boundaryField()[patchi][patchFacei]/146); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.boundaryField()[patchi][patchFacei]*200*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for concrete" [s]

    K_pt.boundaryField()[patchi][patchFacei] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.boundaryField()[patchi][patchFacei],2)) ) * (rho_l*L_v - pc.boundaryField()[patchi][patchFacei]);
}

//- Correct the buildingMaterial lambda (cell)
void Foam::buildingMaterialModels::HamstadConcrete::update_lambda_cell(const volScalarField& w, volScalarField& lambda, label& celli)
{

    lambda.internalField()[celli] = 1.5e0+1.58e-2*w.internalField()[celli];
}

//- Correct the buildingMaterial lambda (boundary)
void Foam::buildingMaterialModels::HamstadConcrete::update_lambda_boundary(const volScalarField& w, volScalarField& lambda, label patchi, label patchFacei)
{

    lambda.boundaryField()[patchi][patchFacei] = 1.5e0+1.58e-2*w.boundaryField()[patchi][patchFacei];
}

//*********************************************************** //
