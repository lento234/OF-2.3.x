// Solves specific humidity transport
{
    DYvt = turbulence->nut()/Sct;
    DYvt.correctBoundaryConditions();

    volScalarField DYvEff("DYvEff", turbulence->nu()/Sc + DYvt);

    fvScalarMatrix YvEqn
    (
        fvm::div(phi, Yv)
      - fvm::laplacian(DYvEff, Yv)
      ==
        vegetation.SYv() // vegetation energy source term
    );

    YvEqn.relax();

    YvEqn.solve();

}
