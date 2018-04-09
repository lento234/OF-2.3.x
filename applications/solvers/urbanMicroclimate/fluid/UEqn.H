    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turb.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
      + vegetation.Su(rho, U) // vegetation momentum source
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve
    (
        UEqn()
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
