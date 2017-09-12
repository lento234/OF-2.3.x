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
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    wnbrName_("undefined-wnbr"),
    TnbrName_("undefined-Tnbr"),
    QrNbrName_("undefined-QrNbr"),
    QrName_("undefined-Qr"),
    QsNbrName_("undefined-QsNbr"),
    QsName_("undefined-Qs")
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
    temperatureCoupledBase(patch(), psf),  
    wnbrName_(psf.wnbrName_),
    TnbrName_(psf.TnbrName_),
    QrNbrName_(psf.QrNbrName_),
    QrName_(psf.QrName_),
    QsNbrName_(psf.QsNbrName_),
    QsName_(psf.QsName_)
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
    temperatureCoupledBase(patch(), dict),
    wnbrName_(dict.lookupOrDefault<word>("wnbr", "none")),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "none")),
    QrNbrName_(dict.lookupOrDefault<word>("QrNbr", "none")),
    QrName_(dict.lookupOrDefault<word>("Qr", "none")),
    QsNbrName_(dict.lookupOrDefault<word>("QsNbr", "none")),
    QsName_(dict.lookupOrDefault<word>("Qs", "none"))    
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
    temperatureCoupledBase(patch(), psf),
    TnbrName_(psf.TnbrName_),
    QrNbrName_(psf.QrNbrName_),
    QrName_(psf.QrName_),
    QsNbrName_(psf.QsNbrName_),
    QsName_(psf.QsName_)	
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

    scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

    const mixedFvPatchScalarField& //CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
            );

    const mixedFvPatchScalarField&
        nbrFieldw = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("w")
            );

    const mixedFvPatchScalarField&
        fieldpc = refCast
            <const mixedFvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("pc")
            );     

    scalarField pcc(fieldpc.patchInternalField()); 

    const fixedValueFvPatchScalarField&
        nbrFieldalphat = refCast
            <const fixedValueFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("alphat")
            );

    const fixedValueFvPatchScalarField&
        nbrFieldmut = refCast
            <const fixedValueFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("mut")
            );            

    const fixedValueFvPatchScalarField&
        nbrFieldgcr = refCast
            <const fixedValueFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("gcr")
            );
    scalarField gcrNbr(nbrFieldgcr.snGrad()/nbrPatch.deltaCoeffs());
    mpp.distribute(gcrNbr);      

////obtain Tambient - can find a better way to import this value?
    Time& time = const_cast<Time&>(nbrMesh.time());
//    label timestep = ceil( (time.value()/3600)-1E-6 ); timestep = timestep%24;

    interpolationTable<scalar> Tambient
    (
	"$FOAM_CASE/0/air/Tambient"
    ); 
	
    interpolationTable<scalar> wambient
    (
	"$FOAM_CASE/0/air/wambient"
    ); 	
///////////

    scalar cap_v = 1880;
    scalar Tref = 273.15; 
    scalar L_v = 2.5e6;                 
    scalar cap_l = 4182;

    // Swap to obtain full local values of neighbour internal field
    scalarField TcNbr(nbrField.patchInternalField()); 
        mpp.distribute(TcNbr);   

    scalarField wcNbr(nbrFieldw.patchInternalField());
    scalar rhoair = 1.2;
    scalarField pv_o = wcNbr*1e5/(0.621945*rhoair);
        mpp.distribute(wcNbr);    
        mpp.distribute(pv_o);               

    scalarField K_pt(Tp.size(), 0.0);
        K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt"); 
        K_pt = (cap_v*(Tc-Tref)+L_v)*K_pt; 

    scalarField lambda_m(Tp.size(), 0.0);
        lambda_m = patch().lookupPatchField<volScalarField, scalar>("lambda_m");                               

    scalarField deltaCoeff_ = nbrPatch.deltaCoeffs(); 
	scalarField mutNbrPatch(nbrFieldmut.patchInternalField() + nbrFieldmut.snGrad()/deltaCoeff_ ); 	
	scalarField alphatNbrPatch(nbrFieldalphat.patchInternalField() + nbrFieldalphat.snGrad()/deltaCoeff_ ); 
	mpp.distribute(mutNbrPatch);  
	mpp.distribute(alphatNbrPatch);  
	mpp.distribute(deltaCoeff_);
    
    scalar cp = 1005; //specific heat of air [J/(kg K)]
    scalar muair = 1.8e-5; scalar Pr = 0.7;
    scalarField heatFlux = (muair/Pr + alphatNbrPatch)*cp*(TcNbr-Tc)*deltaCoeff_; 
//    Info << "heatFlux: " << gSum(heatFlux*patch().magSf()) << endl;
            
    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Tc-5.976*log(Tc));
    scalarField pv_s = pvsat_s*exp((pcc)/(rhol*Rv*Tc));
    
	scalar Dm = 2.5e-5; scalar Sct = 0.7;

    scalarField vaporFlux = (rhoair*Dm + mutNbrPatch/Sct) * (wcNbr-(0.62198*pv_s/1e5)) *deltaCoeff_;         
    scalarField LE = (cap_v*(Tc-Tref)+L_v)*vaporFlux;//Latent and sensible heat transfer due to vapor exchange   */

//    scalarField smoothstep=1/(1+exp((pcc+1000)/30));

    scalarField K_v(Tp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");  
    scalarField Krel(Tp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel");   

    scalarField gl = ((gcrNbr*rhol)/(3600*1000));

//    scalarField CR = gl*cap_l*(Tambient(3600*(timestep+1)) -Tref);
    
	// Calculate rain temperature - approximation for wet-bulb temp///////////
	scalar saturationPressure = 133.322*pow(10,(8.07131-(1730.63/(233.426+Tambient(time.value())))));
	scalar airVaporPressure = wambient(time.value())*1e5/0.621945;
	scalar relhum = airVaporPressure/saturationPressure*100;
	scalar dewPointTemp = Tambient(time.value()) - (100-relhum)/5;
	scalar rainTemp = Tambient(time.value()) - (Tambient(time.value())-dewPointTemp)/3;
	//////////////////////////////////////////////////////////////////////////
	scalarField CR = gl*cap_l*(rainTemp -Tref);

    scalarField Qr(Tp.size(), 0.0);
    if (QrName_ != "none")
    {
        Qr = patch().lookupPatchField<volScalarField, scalar>(QrName_);
    }

    scalarField QrNbr(Tp.size(), 0.0);
    if (QrNbrName_ != "none")
    {
        QrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(QrNbrName_);
        mpp.distribute(QrNbr);
    }
    
    scalarField Qs(Tp.size(), 0.0);
    if (QsName_ != "none")
    {
        Qs = patch().lookupPatchField<volScalarField, scalar>(QsName_);
    }

    scalarField QsNbr(Tp.size(), 0.0);
    if (QsNbrName_ != "none")
    {
        QsNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(QsNbrName_);
        mpp.distribute(QsNbr);
    }   

// term with capillary moisture gradient:                          
    scalarField X = ((cap_l*(Tc-Tref)*Krel)+(cap_v*(Tc-Tref)+L_v)*K_v)*fieldpc.snGrad();
//////////////////////////////////      

    valueFraction() = pos(fieldpc.patchInternalField()+1E3); 
    refValue() = rainTemp;
    refGrad() = (heatFlux + LE + Qr + QrNbr + Qs + QsNbr+CR -X)/(lambda_m+K_pt);
//        Info << "sum(heatFlux): " << sum(heatFlux) << endl;
//        Info << "sum(LE): " << sum(LE) << endl;
//        Info << "sum(CR): " << sum(CR) << endl;
//        Info << "sum(refGrad()):" << sum(refGrad()) << endl;

    mixedFvPatchScalarField::updateCoeffs(); 

    if (debug)
    {
        scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(Tp)
            << " max:" << gMax(Tp)
            << " avg:" << gAverage(Tp)
            << endl;
    } 

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("wnbr")<< wnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("QrNbr")<< QrNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qr")<< QrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("QsNbr")<< QsNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qs")<< QsName_ << token::END_STATEMENT << nl;	
    temperatureCoupledBase::write(os);
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
