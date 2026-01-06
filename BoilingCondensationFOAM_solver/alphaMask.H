Info << "\n--> Executing alphaMask.H with conditional internal update:\n"
     << "    - alpha1 := 0 where alphaMask == 1 (internal) **only if boundary satisfies (alphaMask == 1 && T >= Tact)**\n"
     << "    - alpha1 := 0 where alphaMask == 1 && T >= Tact (boundary)\n"
     << "    at time = " << runTime.timeName() << nl << endl;

// Access required fields
const volScalarField& alphaMask = mesh.lookupObject<volScalarField>("alphaMask");
const volScalarField& T = mesh.lookupObject<volScalarField>("T");
const dimensionedScalar& Tact = mixture->Tact();

// --- BOUNDARY FIELD CHECK ---
bool boundaryConditionSatisfied = false;

forAll(alpha1.boundaryField(), patchI)
{
    const scalarField& maskPatch = alphaMask.boundaryField()[patchI];
    const scalarField& tPatch = T.boundaryField()[patchI];

    forAll(maskPatch, faceI)
    {
        if (maskPatch[faceI] == 1.0 && tPatch[faceI] >= Tact.value())
        {
            boundaryConditionSatisfied = true;
            break;
        }
    }

    if (boundaryConditionSatisfied)
        break;
}

// --- INTERNAL FIELD UPDATE (only if boundary condition is satisfied) ---
if (boundaryConditionSatisfied)
{
    forAll(alpha1.internalField(), cellI)
    {
        if (alphaMask[cellI] == 1.0)
        {
            alpha1[cellI] = 0.0;
        }
    }

    Info << " -> Internal field alpha1 updated (alphaMask == 1)\n";
}
else
{
    Info << " -> Internal field alpha1 NOT updated (boundary condition not met)\n";
}

// --- BOUNDARY FIELD UPDATE (always performed where both mask and T condition are satisfied) ---
forAll(alpha1.boundaryField(), patchI)
{
    const scalarField& maskPatch = alphaMask.boundaryField()[patchI];
    const scalarField& tPatch = T.boundaryField()[patchI];
    scalarField& alpha1Patch = alpha1.boundaryFieldRef()[patchI];

    forAll(alpha1Patch, faceI)
    {
        if (maskPatch[faceI] == 1.0 && tPatch[faceI] >= Tact.value())
        {
            alpha1Patch[faceI] = 0.0;
        }
    }
}

// Ensure alpha2 remains consistent across internal and boundary fields
alpha2.internalFieldRef() = 1.0 - alpha1.internalField();

forAll(alpha2.boundaryField(), patchI)
{
    alpha2.boundaryFieldRef()[patchI] = 1.0 - alpha1.boundaryField()[patchI];
}

Info << " -> alpha1 selectively overwritten (internal: conditional, boundary: based on mask & Tact)\n" << endl;

