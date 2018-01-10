update_sequence = 1; 
pc.dimensions().reset(dimless);
Krel.dimensions().reset(dimless);
K_v.dimensions().reset(dimless);
K_pt.dimensions().reset(dimless);
Crel.dimensions().reset(dimless);
ws.dimensions().reset(dimless);
#include "updatebuildingMaterials.H" 
ws.dimensions().reset(dimensionSet( 1, -3, 0, 0, 0, 0, 0));
Krel.dimensions().reset(dimensionSet( 0, 0, 1, 0, 0, 0, 0));	
K_v.dimensions().reset(dimensionSet( 0, 0, 1, 0, 0, 0, 0));			
K_pt.dimensions().reset(dimensionSet( 1, -1, -1, -1, 0, 0, 0));			
pc.dimensions().reset(dimensionSet( 1, -1, -2, 0, 0, 0, 0));    
Crel.dimensions().reset(dimensionSet( 0, -2, 2, 0, 0, 0, 0)); 

//pc.correctBoundaryConditions();
ws.correctBoundaryConditions();
Krel.correctBoundaryConditions();
Crel.correctBoundaryConditions();
K_v.correctBoundaryConditions();
K_pt.correctBoundaryConditions();//pc.write(); Krel.write(); K_v.write(); w.write(); 

fvScalarMatrix pcEqn
(				
	Crel*fvm::ddt(pc)
	==
	fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
	+fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
);
pcEqn.relax();
pcEqn.solve(); 
