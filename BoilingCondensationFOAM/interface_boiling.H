#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "vector.H"

const scalar deltaAlphaThreshold = 0.4;      // For Δalpha tagging
const scalar cosThetaLimit = 0.866;          // cos(30°)

Info << "\n--> ΔAlpha-Based Interface Detection + Halo Tagging at time = "
     << runTime.timeName() << nl << endl;

// Mixture and fields
autoPtr<temperaturePhaseChangeTwoPhaseMixture> mixture =
    temperaturePhaseChangeTwoPhaseMixture::New(thermo(), mesh);

volScalarField& alpha1 = thermo->alpha1();
volScalarField& int_b = const_cast<volScalarField&>(
    mesh.lookupObject<volScalarField>("int_b")
);

// Create halo marker field
volScalarField haloMarker
(
    IOobject("haloMarker", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Reset fields
// Reset int_b at all patches except "bottom_nucleation"
forAll(int_b.boundaryField(), patchI)
{
    const word& patchName = int_b.boundaryField()[patchI].patch().name();

    if (patchName != "bottom_nucleation" && patchName != "bottom")
    {
        int_b.boundaryFieldRef()[patchI] = 0.0;
    }
}


int_b.internalFieldRef() = 0.0;
haloMarker = 0.0;

// Gradients and mesh topology
volVectorField gradAlpha1 = fvc::grad(alpha1);
volScalarField gradAlphaMag("gradAlphaMag", mag(gradAlpha1));
volVectorField interfaceNormal("interfaceNormal", gradAlpha1 / (gradAlphaMag + dimensionedScalar("tiny", gradAlpha1.dimensions(), VSMALL)));
const labelListList& cellNbrs = mesh.cellCells();
const vectorField& C = mesh.C();

// --- Phase 1: Δalpha-based interface tagging ---
label countInterface = 0;
label countDeltaAlphaTagged = 0;

forAll(alpha1, cellI)
{
    scalar aC = alpha1[cellI];
    const labelList& nbrs = cellNbrs[cellI];
    forAll(nbrs, i)
    {
        scalar aN = alpha1[nbrs[i]];
        if (mag(aC - aN) > deltaAlphaThreshold)
        {
            if (haloMarker[cellI] == 0.0)
            {
                int_b[cellI] = 1.0;
                haloMarker[cellI] = 1.0;
                countInterface++;
                countDeltaAlphaTagged++;
            }
            break;
        }
    }
}
reduce(countInterface, sumOp<label>());
reduce(countDeltaAlphaTagged, sumOp<label>());
Info << "✅ ΔAlpha-based tagging used: " << countDeltaAlphaTagged << " cells\n";

// --- Phase 2: Halo generation (1st and 2nd layer) ---
label countHalo1 = 0, countHalo2 = 0;
forAll(alpha1, cellI)
{
    if (haloMarker[cellI] == 1.0)
    {
        const vector& n = interfaceNormal[cellI];
        const vector& C0 = C[cellI];
        const labelList& nbrs = cellNbrs[cellI];
        scalarField alignment(nbrs.size(), -GREAT);

        forAll(nbrs, i)
        {
            label nbrI = nbrs[i];
            if (haloMarker[nbrI] == 0.0)
            {
                vector d = C[nbrI] - C0;
                scalar cosTheta = (n & d) / (mag(n)*mag(d) + SMALL);
                if (cosTheta > cosThetaLimit)
                    alignment[i] = cosTheta;
            }
        }

        label firstMaxIndex = -1, secondMaxIndex = -1;
        scalar firstMax = -GREAT, secondMax = -GREAT;
        forAll(alignment, i)
        {
            if (alignment[i] > firstMax)
            {
                secondMax = firstMax;
                secondMaxIndex = firstMaxIndex;
                firstMax = alignment[i];
                firstMaxIndex = i;
            }
            else if (alignment[i] > secondMax)
            {
                secondMax = alignment[i];
                secondMaxIndex = i;
            }
        }

        if (firstMaxIndex != -1)
        {
            label h1 = nbrs[firstMaxIndex];
            int_b[h1] = 1.0;
            haloMarker[h1] = 2.0;
            countHalo1++;
        }

        if (secondMaxIndex != -1)
        {
            label h2 = nbrs[secondMaxIndex];
            int_b[h2] = 1.0;
            haloMarker[h2] = 3.0;
            countHalo2++;
        }
    }
}
reduce(countHalo1, sumOp<label>());
reduce(countHalo2, sumOp<label>());


/// --- Phase 3: Boundary tagging with tangential halo ---
label countBoundaryInterface = 0;
label countBoundaryHalo1 = 0;
label countBoundaryHalo2 = 0;

forAll(alpha1.boundaryField(), patchI)
{
    const scalarField& alphaPatch = alpha1.boundaryField()[patchI];
    fvPatchScalarField& int_bPatch = int_b.boundaryFieldRef()[patchI];
    scalarField& haloPatch = haloMarker.boundaryFieldRef()[patchI];

    const label startFace = mesh.boundaryMesh()[patchI].start();
    const vectorField& faceAreas = mesh.boundaryMesh()[patchI].faceAreas();

    forAll(int_bPatch, faceI)
    {
        scalar aFace = alphaPatch[faceI];
        label globalFace = startFace + faceI;
        label owner = mesh.owner()[globalFace];
        scalar aOwner = alpha1[owner];

        if (mag(aOwner - aFace) > deltaAlphaThreshold)
        {
            int_bPatch[faceI] = 1.0;
            haloPatch[faceI] = 1.0;
            countBoundaryInterface++;

            if (haloMarker[owner] == 0.0)
            {
                int_b[owner] = 1.0;
                haloMarker[owner] = 2.0;
                countHalo1++;
            }

            // --- Tangential halo tagging (from owner) ---
            const vector& n = interfaceNormal[owner];
            const vector wallNormal = faceAreas[faceI] / (mag(faceAreas[faceI]) + VSMALL);
            const vector tWall = n - (n & wallNormal) * wallNormal;

            const vector& C0 = C[owner];
            const labelList& nbrs = cellNbrs[owner];

            scalarField alignment(nbrs.size(), -GREAT);
            forAll(nbrs, i)
            {
                label nbr = nbrs[i];
                if (haloMarker[nbr] == 0.0)
                {
                    vector d = C[nbr] - C0;
                    scalar cosTheta = (tWall & d) / (mag(tWall)*mag(d) + SMALL);
                    if (cosTheta > cosThetaLimit)
                        alignment[i] = cosTheta;
                }
            }

            label h1 = -1, h2 = -1;
            scalar max1 = -GREAT, max2 = -GREAT;
            forAll(alignment, i)
            {
                if (alignment[i] > max1)
                {
                    max2 = max1; h2 = h1;
                    max1 = alignment[i]; h1 = nbrs[i];
                }
                else if (alignment[i] > max2)
                {
                    max2 = alignment[i]; h2 = nbrs[i];
                }
            }

            if (h1 != -1)
            {
                int_b[h1] = 1.0;
                haloMarker[h1] = 3.0;
                countBoundaryHalo1++;
            }
            if (h2 != -1)
            {
                int_b[h2] = 1.0;
                haloMarker[h2] = 4.0;
                countBoundaryHalo2++;
            }
        }
    }
}

reduce(countBoundaryInterface, sumOp<label>());
reduce(countBoundaryHalo1, sumOp<label>());
reduce(countBoundaryHalo2, sumOp<label>());

int_b.correctBoundaryConditions();

// --- Final Summary --- //
Info << "\n--- ΔAlpha Interface + Halo Summary ---\n";
Info << " -> Interface method used    : ΔAlpha\n";
Info << " -> Interface cells tagged   : " << countInterface << nl;
Info << " -> 1st-layer halo cells     : " << countHalo1 << nl;
Info << " -> 2nd-layer halo cells     : " << countHalo2 << nl;
Info << " -> Boundary interface faces : " << countBoundaryInterface << nl;
Info << " -> 1st-layer halo at boundary : " << countBoundaryHalo1 << nl;
Info << " -> 2nd-layer halo at boundary : " << countBoundaryHalo2 << nl;

Info << "----------------------------------------\n" << endl;

