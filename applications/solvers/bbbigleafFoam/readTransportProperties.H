    singlePhaseTransportModel laminarTransport(U, phi);

    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta(laminarTransport.lookup("beta"));

    // Reference temperature [K]
    dimensionedScalar TRef(laminarTransport.lookup("TRef"));

    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

    // Laminar Schmidt number
    dimensionedScalar Sc(laminarTransport.lookup("Sc"));

    // Turbulent Schmidt number
    dimensionedScalar Sct(laminarTransport.lookup("Sct"));
