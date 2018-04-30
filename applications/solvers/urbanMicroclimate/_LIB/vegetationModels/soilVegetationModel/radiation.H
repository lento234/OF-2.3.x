
namespace Foam
{

void soilVegetationModel::calc_radiation()
{
    // ---------------------------------------------
    // Direct+Diffused long-wave and short-wave
    // ---------------------------------------------
    
    // Direct+Diffused long-wave and short-wave radiative heat flux at vegetation boundary
    const fvMesh& vegiMesh = mesh_.time().lookupObject<fvMesh>("vegetation"); // "vegetation" a.k.a radiation domain mesh
    const label patchi = vegiMesh.boundaryMesh().findPatchID("air_to_vegetation"); // patch id
    const fvPatch& vegiPatch = vegiMesh.boundary()[patchi]; // patch boundary of vegetation
    scalarField vegiPatchQr = vegiPatch.lookupPatchField<volScalarField, scalar>("Qr"); // long-wave (diffused)
    scalarField vegiPatchQs = vegiPatch.lookupPatchField<volScalarField, scalar>("Qs"); // short-wave (diffused + direct)

    // Net radiation inside vegetation : surface integral
    scalar integrateQr = gSum(vegiPatch.magSf() * vegiPatchQr);
    scalar integrateQs = gSum(vegiPatch.magSf() * vegiPatchQs);
    scalar vegiVolume  = gSum(pos(Cf_.internalField() - 10*SMALL)*mesh_.V());

    Info << "Vegetation : [Radiation]   :: int. Qr = " << integrateQr
         << ", int. Qs = " << integrateQs 
         << ", vol. = " << vegiVolume << endl; 

    // if ((integrateQs/vegiVolume) > (100 * SMALL))
    // {
    //     isDayTime_ = true;
    //     Info << ", it is day time.";
    // }
    // else
    // {
    //     isDayTime_ = false;
    //     Info << ", it is night time.";
    // }
    // Info << endl;
    
    // ---------------------------------------------
    // Direct short-wave radiation through beers law
    // ---------------------------------------------

    // Time of the day
    label timestepsInADay_ = divqrsw.size();
    Time& time = const_cast<Time&>(mesh_.time());
    //Info << "Vegetation : [Radiation]   :: time.value() = " << time.value();
    label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
    //Info << ", 1 timestep = " << timestep;
    timestep = timestep % timestepsInADay_;
    //Info << ", 2 timestep = " << timestep << endl;

    // Net radiation absorbed by vegetation
    scalarList divqrswi =  divqrsw[timestep];

    // radiation density inside vegetation
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > minThreshold)
        {
            qrsw_leaf_[cellI] = (- divqrswi[cellI] )/ LAD_[cellI];
            qrlw_leaf_[cellI] = (integrateQr / vegiVolume)/ LAD_[cellI];
            
        }
        
    }
    qrsw_leaf_.correctBoundaryConditions();
    qrlw_leaf_.correctBoundaryConditions();

    // Integrate direct short-wave radiation absorbed by vegetation
    dimensionedScalar integrateQrsw = fvc::domainIntegrate(qrsw_leaf_ * LAD_);

    //Info << "Vegetation : [Radiation]   :: int. Qrsw = " << integrateQrsw << endl;

    // Missing short-wave radiation, diffused and otherwise.
    scalar missingQs = integrateQs - integrateQrsw.value(); // needs to be added to satisfy energy balance, sadly

    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > minThreshold)
        {
            // Add the missing short-wave radiation
            qrsw_leaf_[cellI] += (missingQs/vegiVolume)/ LAD_[cellI];
            // Calculate the net radiation
            qrad_leaf_[cellI] = qrsw_leaf_[cellI] + qrlw_leaf_[cellI];
            
        }
    }
    qrsw_leaf_.correctBoundaryConditions();
    qrad_leaf_.correctBoundaryConditions();

    Info << "Vegetation : [Radiation]   :: update. int. Qrsw = " << fvc::domainIntegrate(qrad_leaf_ * LAD_) 
         << ", Qr + Qs = " << integrateQr + integrateQs << endl;

}

}
