    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)

      // Does not contain the porous source terms. Porous terms are explicitly added later.
      + turbulence->divDevReff(U)
      //+ fvc::Sp(0.5*Cdf*a*mag(U),U)
      ==
      fvOptions(U)
      //- fvm::Sp(mag(fvc::Sp(Cdf*a, U)), U)
      - fvm::SuSp(0.5*Cdf*a*mag(U), U)

    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
