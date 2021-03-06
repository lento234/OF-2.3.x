/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::radiation::viewFactorUpdated

Description
    View factor radiation model. The system solved is: C q = b
    where:
            Cij  = deltaij/Ej - (1/Ej - 1)Fij
            q    = heat flux
            b    = A eb - Ho
    and:
            eb   = sigma*T^4
            Ej   = emissivity
            Aij  = deltaij - Fij
            Fij  = view factor matrix


SourceFiles
    viewFactorUpdated.C

\*---------------------------------------------------------------------------*/

#ifndef solarLoadModeldirectAndDiffuse_H
#define solarLoadModeldirectAndDiffuse_H

#include "solarLoadModel.H"
#include "singleCellFvMesh.H"
#include "scalarMatrices.H"
#include "globalIndex.H"
#include "scalarListIOList.H"
#include "scalarIOList.H"
#include "mapDistribute.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solarLoad
{

/*---------------------------------------------------------------------------*\
                           Class directAndDiffuse Declaration
\*---------------------------------------------------------------------------*/

class directAndDiffuse
:
    public solarLoadModel
{
    // Private data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Map distributed
        autoPtr<mapDistribute> map_;

        //- Coarse mesh
        singleCellFvMesh coarseMesh_;

        //- Net radiative heat flux [W/m2]
        volScalarField Qs_;

        //- View factor matrix
        autoPtr<scalarSquareMatrix> Fmatrix_;		

        //- Inverse of C matrix
        autoPtr<scalarSquareMatrix> CLU_;
		
		//- skyViewCoeffMatrix_
		autoPtr<scalarRectangularMatrix> skyViewCoeffMatrix_;

        label skyViewCoeffSize; 
		
		//- sunViewCoeffMatrix_
		autoPtr<scalarRectangularMatrix> sunViewCoeffMatrix_;	

        label sunViewCoeffSize;	

        //- Selected patches
        labelList selectedPatches_;

        //- Selected patches
        labelList wallPatchOrNot_;		

        //- Total global coarse faces
        label totalNCoarseFaces_;

        //- Total local coarse faces
        label nLocalCoarseFaces_;
		
        //- Total local wall coarse faces
        label nLocalWallCoarseFaces_;		

        //- Constant emissivity
        bool constEmissivity_;

        //- timestepsInADay
        label timestepsInADay_;

        //- Iterations Counter
        label iterCounter_;

        //- Pivot Indices for LU decomposition
        labelList pivotIndices_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Insert view factors into main matrix
        void insertMatrixElements
        (
            const globalIndex& index,
            const label fromProcI,
            const labelListList& globalFaceFaces,
			const scalarListList& viewFactors,
            scalarSquareMatrix& matrix
        );

        void insertRectangularMatrixElements
        (
            const globalIndex& index,
            const label fromProcI,
            const labelListList& sunskyMap,
            const labelListList& globalFaceFaces,
            const scalarListList& skysunViewCoeffs,
            scalarRectangularMatrix& skysunViewCoeffMatrix
        );          	

        //- Disallow default bitwise copy construct
        directAndDiffuse(const directAndDiffuse&);

        //- Disallow default bitwise assignment
        void operator=(const directAndDiffuse&);


public:

    //- Runtime type information
    TypeName("directAndDiffuse");


    // Constructors

        //- Construct from components
        directAndDiffuse(const volScalarField& T);

        //- Construct from components
        directAndDiffuse(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~directAndDiffuse();


    // Member functions

        // Edit

            //- Solve system of equation(s)
            void calculate();

            //- Read solarLoad properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;


    // Access

        //- Const access to total radiative heat flux field
        inline const volScalarField& Qs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "directAndDiffuseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solarLoad
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
