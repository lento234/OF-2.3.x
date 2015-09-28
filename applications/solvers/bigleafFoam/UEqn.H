    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        fvOptions(U)
      + vegetation.Su(U) // vegetation momentum source
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    // TODO: add buoyancy term later
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
