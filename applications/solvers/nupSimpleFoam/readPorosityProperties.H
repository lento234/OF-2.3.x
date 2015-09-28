    Info<< "Reading porosityProperties\n" << endl;

    IOdictionary porosityProperties
    (
        IOobject
        (
            "porosityProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Assign parameters: Drag coefficient
    dimensionedScalar Cdf
    (
        porosityProperties.lookup("Cdf")
    );
