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

#include "CFDHAMsolidMoistureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
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


CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidMoistureCoupledMixedFvPatchScalarField& psf,
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


CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
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
            "CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::"
            "CFDHAMsolidMoistureCoupledMixedFvPatchScalarField\n"
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


CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidMoistureCoupledMixedFvPatchScalarField& psf,
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

void CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::updateCoeffs()
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

    scalarField pcc(patchInternalField());
    scalarField& pcp = *this;

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);

    const mixedFvPatchScalarField&
        nbrFieldw = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("w")
            );  

/*    const mixedFvPatchScalarField&
        nbrFieldT = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("T")
            ); 
    scalarField TNbr(nbrFieldT.patchInternalField());
    mpp.distribute(TNbr);                         */

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

    const mixedFvPatchScalarField&
        fieldTs = refCast
            <const mixedFvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("Ts")
            );     

    scalarField Ts(fieldTs.patchInternalField());     

    scalarField wcNbr(nbrFieldw.patchInternalField());
    scalar rhoair = 1.2;
//    scalarField pv_o = wcNbr*1e5/(0.621945*rhoair);
        mpp.distribute(wcNbr);
//        mpp.distribute(pv_o);  
//    scalarField pv_o_sat = exp(6.58094e1-7.06627e3/TNbr-5.976*log(TNbr));
//    scalarField pc_o=log(pv_o/pv_o_sat)*rhol*Rv*TNbr; 

    scalarField Krel(pcp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel"); 

//    scalarField Ts(pcp.size(), 0.0);
//        Ts = patch().lookupPatchField<volScalarField, scalar>("Ts");        

    scalarField K_v(pcp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");             

    scalarField deltaCoeff_ = nbrPatch.deltaCoeffs(); 
	scalarField mutNbrPatch(nbrFieldmut.patchInternalField() + nbrFieldmut.snGrad()/deltaCoeff_ ); 
	mpp.distribute(mutNbrPatch);  
	mpp.distribute(deltaCoeff_);
    
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Ts-5.976*log(Ts));
    scalarField pv_s = pvsat_s*exp((pcc)/(rhol*Rv*Ts));

//    scalarField smoothstep=1/(1+exp((pcc+1000)/30));
    scalarField gl = ((gcrNbr*rhol)/(3600*1000));

	scalar Dm = 2.5e-5; scalar Sct = 0.7;

    scalarField vaporFlux = (rhoair*Dm + mutNbrPatch/Sct) * (wcNbr-(0.62198*pv_s/1e5)) *deltaCoeff_; 
//    Info << "vaporFlux: " << gSum(vaporFlux*patch().magSf()) << endl;

// term with temperature gradient:
    scalarField K_pt(pcp.size(), 0.0);
    K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt");                 
	scalarField X = K_pt*fieldTs.snGrad();
//////////////////////////////////                

    valueFraction() = pos(patchInternalField()+1E3);
    refValue() = -1E2;  
    refGrad() = (vaporFlux+gl -X)/(Krel+K_v); 

    mixedFvPatchScalarField::updateCoeffs(); 

    if (debug)
    {
        scalar Q = gSum(kappa(pcp)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(pcp)
            << " max:" << gMax(pcp)
            << " avg:" << gAverage(pcp)
            << endl;
    } 

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::write
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
    CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
