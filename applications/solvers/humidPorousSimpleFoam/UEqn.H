    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
      fvOptions(U)
      - fvm::SuSp(0.5*Cdf*a*mag(U), U) // vegetation drag source term
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    // TODO: add buoyancy term later
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
