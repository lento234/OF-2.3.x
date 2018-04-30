// Check if domain has roots
label hasRoots = 0;
forAll(vegetationRootNames, j)
{
	if (solidRegions[i].name() == vegetationRootNames[j])
		hasRoots = 1;
}

if (hasRoots)
{

	if(pcEqnForm == "mixed")
	{
		//using mixed form for moisture equation - see suGWFoam
		dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
		volScalarField pc_ss = rDeltaT*(ws_old - ws_n); //additional source term appearing due to dw/dt 

		pc.oldTime() = pc_n; //Tricking OpenFOAM to use pc value from previous Picard iteration, instead of previous timestep

		fvScalarMatrix pcEqn
		(
			Crel*fvm::ddt(pc)                
			==
			  fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
			+ fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
			+ pc_ss
			+ vegetation.Sws(Krel,pc) // root uptake
		);
		pcEqn.solve(); 

		pc.oldTime() = pc_old; //restoring pc.oldTime()
	}
	else
	{
		fvScalarMatrix pcEqn
		(
			Crel*fvm::ddt(pc)                
			==
			  fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
			+ fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
			+ vegetation.Sws(Krel,pc) // root uptake
		);
		pcEqn.solve(); 
	}

}
else
{

	if(pcEqnForm == "mixed")
	{
		//using mixed form for moisture equation - see suGWFoam
		dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
		volScalarField pc_ss = rDeltaT*(ws_old - ws_n); //additional source term appearing due to dw/dt 

		pc.oldTime() = pc_n; //Tricking OpenFOAM to use pc value from previous Picard iteration, instead of previous timestep
		
		fvScalarMatrix pcEqn
		(
			Crel*fvm::ddt(pc)                
			==
			  fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
			+ fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
			+ pc_ss
		);
		pcEqn.solve(); 

		pc.oldTime() = pc_old; //restoring pc.oldTime()
	}
	else
	{
		fvScalarMatrix pcEqn
		(
			Crel*fvm::ddt(pc)                
			==
			  fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
			+ fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
		);
		pcEqn.solve(); 
	}

}