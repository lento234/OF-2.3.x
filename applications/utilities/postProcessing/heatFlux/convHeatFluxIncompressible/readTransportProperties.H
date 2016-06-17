
    // Heat capacity
    dimensionedScalar Cp(laminarTransport.lookup("Cp"));

    // Fluid density
    dimensionedScalar rho(laminarTransport.lookup("rho"));

    // Reference Temperature
    dimensionedScalar TRef(laminarTransport.lookup("TRef"));
