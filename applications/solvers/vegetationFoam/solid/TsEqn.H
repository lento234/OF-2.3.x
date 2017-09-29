//update material properties, rho_m and cap_m
update_sequence = 2;
#include "updatebuildingMaterials.H"

rho_m.correctBoundaryConditions();
cap_m.correctBoundaryConditions();
lambda_m.correctBoundaryConditions();

// dimensionedScalar cap_l("cap_l",dimensionSet(0,2,-2,-1,0,0,0),scalar(4182));
// dimensionedScalar cap_v("cap_v",dimensionSet(0,2,-2,-1,0,0,0),scalar(1880));
dimensionedScalar Tref("Tref",dimensionSet(0,0,0,1,0,0,0),scalar(273.15));
// dimensionedScalar L_v("L_v",dimensionSet(0,2,-2,0,0,0,0), 2.5e6);//   +(cap_l.value()-cap_v.value())*Tref.value());

//volScalarField C_t = rho_m*cap_m + ws*cap_l;
volScalarField C_t = rho_m*cap_m;

fvScalarMatrix TsEqn
(
	C_t*fvm::ddt(Ts)
	==
	fvm::laplacian(lambda_m,Ts,"laplacian(Krel,pc)")
	// + fvm::laplacian(((Ts-Tref)*cap_v+L_v)*K_pt,Ts,"laplacian(Krel,pc)")
	// + fvc::laplacian((Ts-Tref)*cap_l*Krel,pc,"laplacian(Krel,pc)")
	// + fvc::laplacian(((Ts-Tref)*cap_v+L_v)*K_v,pc,"laplacian(Krel,pc)")
);
TsEqn.relax();
TsEqn.solve();
