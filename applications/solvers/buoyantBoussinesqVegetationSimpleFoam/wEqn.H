// Solves specific humidity transport
{
    Dwt = turbulence->nut()/Sct;
    Dwt.correctBoundaryConditions();

    volScalarField DwEff("DwEff", turbulence->nu()/Sc + Dwt);

    fvScalarMatrix wEqn
    (
        fvm::div(phi, w)
      - fvm::laplacian(DwEff, w)
      ==
        vegetation.Sw() // vegetation energy source term
    );

    wEqn.relax();

    wEqn.solve();

}
