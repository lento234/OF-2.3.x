{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        fvOptions(T)
      +  vegetation.Sh() // vegetation energy source term
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);
}
