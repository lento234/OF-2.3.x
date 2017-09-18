/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::buildingMaterialModels::Soil

Description
    buildingMaterial constitutive relationship as in Soil et al. (1977).
   
    Reference:
    Soil R, Vauclin M, Touma J, Wierenga P, Vachaud G. Comparision of numerical
    simulation models for one dimensional infiltration. buildingMaterial Sci. Soc. Am J. 1977,
    41:285-294

SourceFiles
    Soil.C

\*---------------------------------------------------------------------------*/

#ifndef Soil_H
#define Soil_H

#include "buildingMaterialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{

/*---------------------------------------------------------------------------*\
                           Class Soil Declaration
\*---------------------------------------------------------------------------*/

class Soil
:
    public buildingMaterialModel
{
    // Private data

        dictionary SoilCoeffs_;

        dimensionedScalar rho_;
        dimensionedScalar cap_;    

        //dimensionedTensor Ks_;
        //dimensionedScalar theta_s_;
        //dimensionedScalar theta_r_;
        //dimensionedScalar alpha_;
        //dimensionedScalar beta_;
        //dimensionedScalar gamma_;
        //dimensionedScalar A_;
        //dimensionedScalar Ss_; 

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("Soil");


    // Constructors

        //- Construct from components
        Soil
        (
            const word& name, 
            const dictionary& buildingMaterialProperties,
            const word& cellZoneModel
            //volScalarField& h,
            //volScalarField& theta,
            //volScalarField& kr,
            //volScalarField& Ch
        );


    // Destructor

        ~Soil()
        {}


    // Member Functions

        //- Correct the buildingMaterial moisture content (cell)
        void update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli);

        //- Correct the buildingMaterial moisture content (boundary)
        void update_w_C_boundary(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label patchi, label patchFacei); 

        //- Correct the buildingMaterial liquid permeability (cell)
        void update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli);

        //- Correct the buildingMaterial liquid permeability (boundary)
        void update_Krel_boundary(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label patchi, label patchFacei); 

        //- Correct the buildingMaterial vapor permeability (cell)
        void update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli);

        //- Correct the buildingMaterial vapor permeability (boundary)
        void update_Kv_boundary(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label patchi, label patchFacei); 

        //- Correct the buildingMaterial K_pt (cell)
        void update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli);

        //- Correct the buildingMaterial K_pt (boundary)
        void update_Kpt_boundary(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label patchi, label patchFacei);

        //- Correct the buildingMaterial lambda (cell)
        void update_lambda_cell(const volScalarField& w, volScalarField& lambda, label& celli);

        //- Correct the buildingMaterial lambda (boundary)
        void update_lambda_boundary(const volScalarField& w, volScalarField& lambda, label patchi, label patchFacei);                      

        //- Read transportProperties dictionary
        bool read(const dictionary& buildingMaterialProperties);

        dimensionedScalar rho()
        {
            return rho_;
        }

        dimensionedScalar cap()
        {
            return cap_;
        }    

        /*dimensionedTensor Ks()
        {
            return Ks_;
        }

        dimensionedScalar theta_s()
        {
            return theta_s_;
        }

        dimensionedScalar Ss()
        {
            return Ss_;
        }*/

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace buildingMaterialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
