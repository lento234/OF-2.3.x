/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    QrNbrName_("undefined-QrNbr"),
    QsNbrName_("undefined-QsNbr"),
    impermeable_("undefined-impermeable")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    QrNbrName_(psf.QrNbrName_),
    QsNbrName_(psf.QsNbrName_),
    impermeable_(psf.impermeable_)
{}


CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    QrNbrName_(dict.lookupOrDefault<word>("QrNbr", "none")),
    QsNbrName_(dict.lookupOrDefault<word>("QsNbr", "none")),
    impermeable_(dict.lookupOrDefault<bool>("impermeable", false))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::"
            "CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }  

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    QrNbrName_(psf.QrNbrName_),
    QsNbrName_(psf.QsNbrName_),
    impermeable_(psf.impermeable_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];

    scalar cap_v = 1880; scalar Tref = 273.15; scalar L_v = 2.5e6; scalar cap_l = 4182;
    scalar cp = 1005; //specific heat of air [J/(kg K)]
    scalar muair = 1.8e-5; scalar Pr = 0.7;
    scalar Dm = 2.5e-5; scalar Sct = 0.7;
    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);

    //scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

    const mixedFvPatchScalarField& //CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("T")
            );
    scalarField TcNbr(nbrField.patchInternalField()); 
        mpp.distribute(TcNbr);
    scalarField TNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("T");
        mpp.distribute(TNbr);

    const mixedFvPatchScalarField&
        nbrFieldw = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("w")
            );
    scalarField wcNbr(nbrFieldw.patchInternalField());
        mpp.distribute(wcNbr);
    scalarField wNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("w");
    scalarField rhoNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("rho");
    scalarField pv_o = wNbr*1e5/(0.621945*rhoNbr);
        mpp.distribute(wNbr);
        mpp.distribute(rhoNbr);  
        mpp.distribute(pv_o); 

    const mixedFvPatchScalarField&
        fieldpc = refCast
            <const mixedFvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("pc")
            );
    scalarField pc(Tp.size(), 0.0);
        pc = patch().lookupPatchField<volScalarField, scalar>("pc");    

    scalarField K_pt(Tp.size(), 0.0);
        K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt"); 
        K_pt = (cap_v*(Tp-Tref)+L_v)*K_pt; 

    scalarField lambda_m(Tp.size(), 0.0);
        lambda_m = patch().lookupPatchField<volScalarField, scalar>("lambda_m");                               

    scalarField deltaCoeff_ = nbrPatch.deltaCoeffs();
        mpp.distribute(deltaCoeff_);
    scalarField alphatNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("alphat");
        mpp.distribute(alphatNbr);
    scalarField mutNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("mut");
        mpp.distribute(mutNbr); 
    
    scalarField heatFlux = (muair/Pr + alphatNbr)*cp*(TcNbr-Tp)*deltaCoeff_; 
            
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Tp-5.976*log(Tp));
    scalarField pv_s = pvsat_s*exp((pc)/(rhol*Rv*Tp));
    
    scalarField vaporFlux = (rhoNbr*Dm + mutNbr/Sct) * (wcNbr-(0.62198*pv_s/1e5)) *deltaCoeff_;         
    scalarField LE = (cap_v*(Tp-Tref)+L_v)*vaporFlux;//Latent and sensible heat transfer due to vapor exchange   */

    scalarField K_v(Tp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");  
    scalarField Krel(Tp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel");   

    scalarField gcrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("gcr");
        mpp.distribute(gcrNbr); 

    scalarField gl = ((gcrNbr*rhol)/(3600*1000));

    // Calculate rain temperature - approximation for wet-bulb temp///////////
    //obtain Tambient - can find a better way to import this value?
    Time& time = const_cast<Time&>(nbrMesh.time());
    //label timestep = ceil( (time.value()/3600)-1E-6 ); timestep = timestep%24;

    interpolationTable<scalar> Tambient
    (
        "$FOAM_CASE/0/air/Tambient"
    ); 
    
    interpolationTable<scalar> wambient
    (
        "$FOAM_CASE/0/air/wambient"
    );     
    ///////////
    scalar saturationPressure = 133.322*pow(10,(8.07131-(1730.63/(233.426+Tambient(time.value())))));
    scalar airVaporPressure = wambient(time.value())*1e5/0.621945;
    scalar relhum = airVaporPressure/saturationPressure*100;
    scalar dewPointTemp = Tambient(time.value()) - (100-relhum)/5;
    scalar rainTemp = Tambient(time.value()) - (Tambient(time.value())-dewPointTemp)/3;
    //////////////////////////////////////////////////////////////////////////

    scalarField QrNbr(Tp.size(), 0.0);
    if (QrNbrName_ != "none")
    {
        QrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(QrNbrName_);
        mpp.distribute(QrNbr);
    }

    scalarField QsNbr(Tp.size(), 0.0);
    if (QsNbrName_ != "none")
    {
        QsNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(QsNbrName_);
        mpp.distribute(QsNbr);
    }   

        //-- Gravity-enthalpy flux --//
        //lookup gravity vector
        uniformDimensionedVectorField g = db().lookupObject<uniformDimensionedVectorField>("g");
        scalarField gn = g.value() & patch().nf();

        scalarField phiGT = (cap_l*(Tp-Tref))*Krel*rhol*gn;


        // term with capillary moisture gradient:                          
        scalarField X = ((cap_l*(Tp-Tref)*Krel)+(cap_v*(Tp-Tref)+L_v)*K_v)*fieldpc.snGrad();
        //////////////////////////////////  
        scalarField CR = ( pos(patchInternalField()+fieldpc.snGrad()/patch().deltaCoeffs()+1E3)*((Krel+K_v)*fieldpc.snGrad() - Krel*rhol*gn) +
                           neg(patchInternalField()+fieldpc.snGrad()/patch().deltaCoeffs()+1E3)*gl ) * cap_l*(rainTemp -Tref) * pos(gl-VSMALL);

        

    if(impermeable_ == false)
    {
        valueFraction() = 0;//pos(fieldpc.patchInternalField()+1E3); 
        refValue() = 0;//rainTemp;
        refGrad() = (heatFlux + LE + QrNbr + QsNbr + CR + phiGT -X)/(lambda_m+K_pt);
//      refGrad() = (heatFlux + LE + QrNbr + QsNbr + CR -X)/(lambda_m+K_pt);
    }
    else
    {
        valueFraction() = 0;
        refValue() = 0;
        refGrad() = (heatFlux + QrNbr + QsNbr)/(lambda_m);
    }
//Info << "111: " << heatFlux << " " << LE << " " << -X << " " << K_pt << " " << refGrad() << endl;
    mixedFvPatchScalarField::updateCoeffs(); 

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("QrNbr")<< QrNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("QsNbr")<< QsNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("impermeable")<< impermeable_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
