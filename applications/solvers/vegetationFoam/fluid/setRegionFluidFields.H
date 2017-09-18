    const fvMesh& mesh = fluidRegions[i];

    rhoThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    radiation::radiationModel& rad = radiation[i];
    solarLoad::solarLoadModel& sol = solarLoad[i];

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& w = wFluid[i];
    //volScalarField& gcr = gcrFluid[i];
    volScalarField& Tl = TlFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    fv::IOoptionList& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    vegetationModel &vegetation = vegetationFluid[i];
