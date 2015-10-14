    Info<< "Reading vegetationProperties\n" << endl;

    IOdictionary vegetationProperties
    (
        IOobject
        (
            "vegetationProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Assign parameters: Drag coefficient
    dimensionedScalar Cdf
    (
        vegetationProperties.lookup("Cdf")
    );
