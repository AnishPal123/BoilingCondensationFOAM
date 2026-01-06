/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "constant.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "twoPhaseMixtureEThermo.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        temperaturePhaseChangeTwoPhaseMixture,
        constant,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::constant
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    temperaturePhaseChangeTwoPhaseMixture(mixture, mesh),
    coeffC_
    (
        "coeffC",
        dimless/dimTime/dimTemperature,
        optionalSubDict(type() + "Coeffs")
    ),
    coeffE_
    (
        "coeffE",
        dimless/dimTime/dimTemperature,
        optionalSubDict(type() + "Coeffs")
    ),
    
      Tact_
    (
    "Tact",
    dimTemperature,
    optionalSubDict(type() + "Coeffs")
   ) ,
   
     betaL_
   (
    "betaL",
    dimless/dimTemperature,
    optionalSubDict(type() + "Coeffs")
   ),
   
    betav_
   (
    "betav",
    dimless/dimTemperature,
    optionalSubDict(type() + "Coeffs")
   )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//returns Tact
const Foam::dimensionedScalar&
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::Tact() const
{
    return Tact_;
}

//returns betaL
const Foam::dimensionedScalar&
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::betaL() const
{
    return betaL_;
}
//returns betav
const Foam::dimensionedScalar&
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::betav() const
{
    return betav_;
}


//returns TSat
const Foam::dimensionedScalar&
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::TSat() const
{
    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();
    return TSat;
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotAlphal() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const volScalarField& int_b = mesh_.lookupObject<volScalarField>("int_b");
    const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");
    

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        Cc_mod*int_b*coeffC_*mixture_.rho2()*max(TSat - T, T0),
       -coeffE_*int_b*mixture_.rho1()*max(T - TSat, T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDot() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
     const volScalarField& int_b = mesh_.lookupObject<volScalarField>("int_b");
    const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    volScalarField mDotE
    (
        "mDotE",
        coeffE_*int_b*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
      * max(T - TSat, T0)
    );
    volScalarField mDotC
    (
        "mDotC",
        Cc_mod*int_b*coeffC_*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
      * max(TSat - T, T0)
    );

    if (mesh_.time().writeTime())
    {
        mDotC.write();
        mDotE.write();
    }

    return Pair<tmp<volScalarField>>
    (
        tmp<volScalarField>(new volScalarField(mDotC)),
        tmp<volScalarField>(new volScalarField(-mDotE))
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotDeltaT() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
     const volScalarField& int_b = mesh_.lookupObject<volScalarField>("int_b");
    const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    return Pair<tmp<volScalarField>>
    (
        (
            coeffC_*Cc_mod*int_b*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
          * pos(TSat - T)
        ),
        (
            coeffE_*int_b*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
          * pos(T - TSat)
        )
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::TSource() const
{

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
     const volScalarField& int_b = mesh_.lookupObject<volScalarField>("int_b");
    const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");

    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime
        )
    );

    fvScalarMatrix& TSource = tTSource.ref();

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar L = mixture_.Hf2() - mixture_.Hf1();

    const volScalarField Vcoeff
    (
        coeffE_*int_b*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
      * L*pos(T - TSat)
    );
    const volScalarField Ccoeff
    (
        Cc_mod*int_b*coeffC_*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
      * L*pos(TSat - T)
    );

    TSource =
        fvm::Sp(Vcoeff, T) - Vcoeff*TSat
      + fvm::Sp(Ccoeff, T) - Ccoeff*TSat;

    return tTSource;
}


void Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::correct()
{
}


bool Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::read()
{
    if (temperaturePhaseChangeTwoPhaseMixture::read())
    {
        subDict(type() + "Coeffs").readEntry("coeffC", coeffC_);
        subDict(type() + "Coeffs").readEntry("coeffE", coeffE_);
        subDict(type() + "Coeffs").readEntry("Tact", Tact_);
         subDict(type() + "Coeffs").readEntry("betaL", betaL_);
        subDict(type() + "Coeffs").readEntry("betav", betav_);

        return true;
    }

    return false;
}


// ************************************************************************* //
