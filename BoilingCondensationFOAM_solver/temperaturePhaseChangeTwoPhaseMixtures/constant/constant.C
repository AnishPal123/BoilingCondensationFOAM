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
    Cc_
    ("Cc", 
    dimless, 
    optionalSubDict(type() + "Coeffs")
    ),
    
   Cv_
   ("Cv", 
   dimless, 
   optionalSubDict(type() + "Coeffs")
   ),
   
    n_
    ("n", 
    dimless/dimVolume, 
    optionalSubDict(type() + "Coeffs")
    ),
    
    dNuc_
    ("dNuc", 
    dimLength, 
    optionalSubDict(type() + "Coeffs")
    ),
    
    PSat_
   (
    "PSat",
    dimPressure,
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
   
 
   

{
  correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::rRb
(
    const volScalarField& limitedAlpha1
) const
{
    return pow
    (
        ((4*Foam::constant::mathematical::pi*n_)/3)
       *limitedAlpha1/(1.0 + alphaNuc() - limitedAlpha1),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*(Foam::constant::mathematical::pi)*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}

Foam::tmp<Foam::volScalarField>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::pCoeff
    (
        const volScalarField& p,
        const volScalarField& T
    ) const
    

{
    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );
    
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );
    

    
    volScalarField rho
    (
        limitedAlpha1*mixture_.rho1() + (scalar(1) - limitedAlpha1)*mixture_.rho2()
    );
    
    return
        (3*mixture_.rho2()*mixture_.rho1())*sqrt(2/(3*mixture_.rho1()))
        *rRb(limitedAlpha1)/(rho*sqrt(mag(p - PSat_) + SMALL*PSat_));
}


//returns PSat
const Foam::dimensionedScalar&
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::PSat() const
{
    return PSat_;
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

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotAlphal() const
{
	const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
        const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
        const volScalarField& int_b = mesh_.lookupObject<volScalarField>("int_b");
        const volScalarField& F_Cav = mesh_.lookupObject<volScalarField>("F_Cav");
        const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");
       
       
               volScalarField pCoeff(this->pCoeff(p, T));

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);
    const dimensionedScalar p0(dimPressure, Zero);
    
     volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );      
    


    return Pair<tmp<volScalarField>>
    (
        Cc_*F_Cav*limitedAlpha1*pCoeff*max(p - PSat_, p0) + Cc_mod*coeffC_*int_b*mixture_.rho2()*max(TSat - T, T0),
       Cv_*F_Cav*(1.0 + alphaNuc()-limitedAlpha1)*pCoeff*min(p - PSat_, p0)-coeffE_*int_b*mixture_.rho1()*max(T - TSat, T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotT() const
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
        coeffE_*int_b*mixture_.rho1()*int_b*clamp(mixture_.alpha1(), zero_one{})
      * max(T - TSat, T0)
    );
    volScalarField mDotC
    (
        "mDotC",
        Cc_mod*coeffC_*int_b*mixture_.rho2()*int_b*clamp(mixture_.alpha2(), zero_one{})
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
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotP() const
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
                const volScalarField& F_Cav = mesh_.lookupObject<volScalarField>("F_Cav");
                    const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");
    
    
    volScalarField pCoeff(this->pCoeff(p, T));

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    /*const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );*/

    //const volScalarField PSat_F = thermo.PSat_F(this->TPSat_F(T));
    volScalarField apCoeff(limitedAlpha1*pCoeff);
    return Pair<tmp<volScalarField>>
    (
       Cc_*F_Cav*(scalar(1) - limitedAlpha1)*pos0(p - PSat_)*apCoeff,
       (-Cv_)*F_Cav*(1.0 + alphaNuc() - limitedAlpha1)*neg(p - PSat_)*apCoeff    
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
            Cc_mod*coeffC_*int_b*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
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
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volScalarField& int_b = mesh_.lookupObject<volScalarField>("int_b");
        const volScalarField& F_Cav = mesh_.lookupObject<volScalarField>("F_Cav");
                    const volScalarField& Cc_mod = mesh_.lookupObject<volScalarField>("Cc_mod");
        
    const dimensionedScalar p0(dimPressure, Zero);

    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime
        )
    );

    fvScalarMatrix& TSource = tTSource.ref();
            volScalarField pCoeff(this->pCoeff(p, T));

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar L = mixture_.Hf2() - mixture_.Hf1();
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    const volScalarField Vcoeff
    (
        coeffE_*int_b*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
      * L*pos(T - TSat)
    );
    const volScalarField Ccoeff
    (
        Cc_mod*coeffC_*int_b*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
      * L*pos(TSat - T)
    );
     const volScalarField VcoeffP
    (
       Cv_*F_Cav*alphaNuc()*limitedAlpha1*pCoeff*min(p - PSat_, p0)*L
    );
   
   Info << "dim VcoeffP: " << (Cv_ * alphaNuc() * limitedAlpha1 * pCoeff * min(p - PSat_, p0) * L).ref().dimensions() << nl;

    const volScalarField CcoeffP
    (
      Cc_*F_Cav*(scalar(1) - limitedAlpha1)*pCoeff*max(p - PSat_, p0)*L
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
        subDict(type() + "Coeffs").readEntry("Cc", Cc_);
        subDict(type() + "Coeffs").readEntry("Cv", Cv_);
        subDict(type() + "Coeffs").readEntry("n", n_);
        subDict(type() + "Coeffs").readEntry("dNuc", dNuc_);
        subDict(type() + "Coeffs").readEntry("PSat", PSat_);
        subDict(type() + "Coeffs").readEntry("Tact", Tact_);
         subDict(type() + "Coeffs").readEntry("betaL", betaL_);
        subDict(type() + "Coeffs").readEntry("betav", betav_);
        
        return true;
    }

    return false;
}


// ************************************************************************* //
