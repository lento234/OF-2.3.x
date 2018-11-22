    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
      fvOptions(U)
      - fvm::SuSp(cd*a*mag(U), U) // should be twice as much
      //- fvm::SuSp(0.5*cd*a*mag(U), U)

    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
