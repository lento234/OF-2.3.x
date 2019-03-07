{
    dimensionedScalar Dm("Dm",dimensionSet(0,2,-1,0,0,0,0),scalar(2.5e-5)); 
    scalar Sct = 0.7;

    fvScalarMatrix wEqn
    (
        fvm::ddt(rho,w)
      + fvm::div(phi, w)
//      - fvm::laplacian(turb.alphaEff(), w) //In reactingParcelFimFoam, they use alphaEff, but in others Mass Diffusivity = muEff/Sct      
      - fvm::laplacian(rho*Dm + turb.mut()/Sct, w)
//      - fvm::laplacian(rho*Dm, w)      
    );
    wEqn.relax();

    wEqn.solve();

}
