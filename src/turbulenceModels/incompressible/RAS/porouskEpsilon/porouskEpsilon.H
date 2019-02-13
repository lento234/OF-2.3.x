/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::porouskEpsilon

Authour
    Lento Manickathan. <manickathan@arch.ethz.ch>

Description
    Standard k-epsilon turbulence model for incompressible flows.

    Reference:
    \verbatim
        Launder, B.E., and Spalding, D.B.,
        "Mathematical Models of Turbulence",
        Academic Press, 1972.

        Launder, B.E., and Spalding, D.B.,
        "Computational Methods for Turbulent Flows",
        Comp. Meth in Appl Mech & Eng'g, Vol 3, 1974, pp 269-289.

        Liang, L., Xiaofeng, L., Borong, L., & Yingxin, Z. (2006).
        "Improved $k-\epsilon$ two-equation turbulence model for canopy flow.
        Atmospheric Environment, 40(4), 762–770.

        Gromke, C., Blocken, B., Janssen, W., Merema, B.,
        van Hooff, T., & Timmermans, H. (2014).
        "CFD analysis of transpirational cooling by vegetation: Case study for
        specific meteorological conditions during a heat wave in Arnhem, Netherlands."
        Building and Environment, 83, 11–26.

    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        porouskEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            sigmaK      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    porouskEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef porouskEpsilon_H
#define porouskEpsilon_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class porouskEpsilon Declaration
\*---------------------------------------------------------------------------*/

class porouskEpsilon
:
    public RASModel
{

protected:

    // Protected data

        // Flag to turn model options on

            // Switch buoyancyProduction_;


            // word TName_;

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar cd_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar betaP_;
            dimensionedScalar betaD_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField a_;
            volScalarField nut_;

public:

    //- Runtime type information
    TypeName("porouskEpsilon");

    // Constructors

        //- Construct from components
        porouskEpsilon
        (
            const volVectorField& U,
            //const volScalarField& LAD,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~porouskEpsilon()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- return LAD
        virtual tmp<volScalarField> a() const
        {
            return a_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
