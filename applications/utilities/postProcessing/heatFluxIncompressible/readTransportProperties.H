
  // Heat capacity
  dimensionedScalar Cp(laminarTransport.lookup("Cp"));

  // Thermal expansion coefficient [1/K]
  dimensionedScalar Pr(laminarTransport.lookup("Pr"));

  // Turbulent Prandtl number
  dimensionedScalar Prt(laminarTransport.lookup("Prt"));

  // Fluid density
  dimensionedScalar rho(laminarTransport.lookup("rho"));

  // Reference Temperature
  dimensionedScalar TRef(laminarTransport.lookup("TRef"));
