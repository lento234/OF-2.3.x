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

    if (transientSimple.momentumPredictor())
    {
        solve
        (
            UEqn()
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
