// Solves specific humidity transport
{
    Dqt = turbulence->nut()/Sct;
    Dqt.correctBoundaryConditions();

    volScalarField DqEff("DEff", turbulence->nu()/Sc + Dqt);

    fvScalarMatrix qEqn
    (
        fvm::div(phi, q)
      - fvm::laplacian(DqEff, q)
      ==
        vegetation.Sq(U, T, q) // vegetation energy source term
    );

    qEqn.relax();

    qEqn.solve();

}
