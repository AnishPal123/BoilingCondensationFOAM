/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::phaseChangeTwoPhaseMixtures::constant

Description
    constant condensation/saturation model.


SourceFiles
    constant.C

\*--------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "temperaturePhaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                             Class constant
\*--------------------------------------------------------------------*/

class constant
:
    public temperaturePhaseChangeTwoPhaseMixture
{
    // Private data

        //- Condensation coefficient [1/s/K]
        dimensionedScalar coeffC_;

        //- Evaporation coefficient [1/s/K]
        dimensionedScalar coeffE_;
        
          //- Condensation rate coefficient [-]
          dimensionedScalar Cc_;

        //- Vapourisation rate coefficient [-]
        dimensionedScalar Cv_;
        
        //- Bubble number density
        dimensionedScalar n_;
        
        

        //- Nucleation site diameter
        dimensionedScalar dNuc_;
        
          //- Sat Pressure (PSat_)
        dimensionedScalar PSat_;
        
        //-Activation Temp (Tact_)
         dimensionedScalar Tact_;
         //-betaL and betav
         dimensionedScalar betaL_;
         dimensionedScalar betav_;
       //- Nucleation site volume-fraction
        dimensionedScalar alphaNuc() const;

        //- Reciprocal bubble radius
        tmp<volScalarField>rRb(const volScalarField& limitedAlpha1) const;
        
       //- Part of the condensation and vaporisation rates [s/m2]
        tmp<volScalarField> pCoeff
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constant() = default;


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const;
        //- Return the pressure term of the mass condensation and vaporisation 
        //  rates as coefficients to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const;

        //- Return the temperature term of the mass condensation and 
        //  vaporisation rates as coefficients
        virtual Pair<tmp<volScalarField>> mDotT() const;
        
        
        //- Return saturation pressure as a dimensioned scalar
       virtual const dimensionedScalar& PSat() const override;
       //- Return activation Tact as a dimensioned scalar
       virtual const dimensionedScalar& Tact() const override;
        virtual const dimensionedScalar& betaL() const override;
       virtual const dimensionedScalar& betav() const override;
       virtual const dimensionedScalar& TSat() const override;
        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (Tsat - T) for the condensation rate
        //  and a coefficient to multiply  (T - Tsat) for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotDeltaT() const;

        //- Source for T equarion
        virtual tmp<fvScalarMatrix> TSource() const;

        //- Correct the constant phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace temperaturePhaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
