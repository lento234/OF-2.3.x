/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::vegetationModel

Description
    Namespace for vegetation modelling

Author
    Lento Manickathan

SourceFiles
    vegetationModel.C
    vegetationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef vegetationModel_H
#define vegetationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvm.H"
#include "bound.H"
#include "Switch.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{
/*---------------------------------------------------------------------------*\
                         Class vegetationModel Declaration
\*---------------------------------------------------------------------------*/

class vegetationModel
:
    public IOdictionary

{
protected:

    // Protected Data

        //- Reference to vegetation model dictionary
        dictionary vegetationDict_;
    
        //- Reference to mesh database
        const fvMesh& mesh_;

        //- Reference to time database
        const Time& runTime_;

        //- Leaf area density
        //volScalarField LAD_;    // leaf area density
        
private:
    // disallow copy construct
    vegetationModel(const vegetationModel&);

    // disallow default bitwise assignment
    void operator=(const vegetationModel&);


public:
    
    //- Runtime type info
    TypeName("vegetationModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            vegetationModel,
            dictionary,
            (
                const dictionary& vegetationDict
                //const volScalarField& LAD
            ),
            (vegetationDict)
        );

    // Selectors

        //- Return a reference to the selected buildingMaterial model
        static autoPtr<vegetationModel> New
        (
            const word& name,
            const dictionary& vegetationDict
        );

    // Constructors

        //- Construct from components
        vegetationModel
        (
            const word& name,
            const dictionary& vegetationDict
        );


     // Destructor

        virtual ~vegetationModel()
        {}


    // Member Functions

        // -----------------------------------------------

        // return vegetation properties dictionary
        dictionary vegetationDict() const
        {
            return vegetationDict_;
        }

        // -----------------------------------------------

        
        // solve all
        virtual void solve();

        // -----------------------------------------------
        // source terms

        // return vegetation energy source
        virtual tmp<volScalarField> Sh() = 0;

        // return vegetation momentum source
        virtual tmp<fvVectorMatrix> Su(volVectorField& U) = 0;

        // return vegetation specific humidity source
        virtual tmp<volScalarField> Sq() = 0;

};


} // end namepsace Foam

#endif

// ************************************************************************* //
