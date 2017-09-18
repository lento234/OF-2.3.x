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

#include "CFDHAMfluidMoistureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
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


CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMfluidMoistureCoupledMixedFvPatchScalarField& psf,
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


CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
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
            "CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::"
            "CFDHAMfluidMoistureCoupledMixedFvPatchScalarField\n"
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


CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMfluidMoistureCoupledMixedFvPatchScalarField& psf,
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

void CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::updateCoeffs()
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

//    scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

/*    const mixedFvPatchScalarField& //CFDHAMfluidMoistureCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(wnbrName_)
            );   */

    const mixedFvPatchScalarField& //CFDHAMfluidMoistureCoupledMixedFvPatchScalarField&
        nbrTField = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
            );   

    const mixedFvPatchScalarField& //CFDHAMfluidMoistureCoupledMixedFvPatchScalarField&
        nbrpcField = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("pc")
            );                       

    // Swap to obtain full local values of neighbour internal field
//    scalarField wcNbr(nbrField.patchInternalField());
//    mpp.distribute(wcNbr);

    scalarField TNbr(nbrTField.patchInternalField());
    mpp.distribute(TNbr);  

    scalarField pcNbr(nbrpcField.patchInternalField());
    mpp.distribute(pcNbr); 

    scalarField p(Tp.size(), 0.0);
        p = patch().lookupPatchField<volScalarField, scalar>("p");    

    scalarField rhoair(Tp.size(), 0.0);
        rhoair = patch().lookupPatchField<volScalarField, scalar>("rho");            

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/TNbr-5.976*log(TNbr));
    
    scalarField pv_s = pvsat_s*exp((pcNbr)/(rhol*Rv*TNbr));
    
    valueFraction() = 1.0;//KDeltaNbr/(KDeltaNbr + KDelta);
    refValue() = 0.62198*pv_s/p;//pv_s/pvsat_s;
    refGrad() = 0.0;//(Qr + QrNbr + Qs + QsNbr)/(kappa(Tp));

    mixedFvPatchScalarField::updateCoeffs();

/*    if (debug)
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
    } */

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::write
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
    CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
