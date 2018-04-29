{
    dimensionedScalar Dc("Dc", dimensionSet(0,2,-1,0,0,0,0), scalar(1.6e-5)); // Diffusivity of CO2
    scalar Sct = 0.7; // equal to vapor (but shouldn't be higher?)

    fvScalarMatrix cEqn
    (
        fvm::ddt(rho, c)
      + fvm::div(phi, c)
      - fvm::laplacian(rho*Dc + turb.mut()/Sct, c)
      ==
        vegetation.Sc() // vegetation CO2 source term (kg/m3 s * mol/mol)
    );
    cEqn.relax();

    cEqn.solve();

}
