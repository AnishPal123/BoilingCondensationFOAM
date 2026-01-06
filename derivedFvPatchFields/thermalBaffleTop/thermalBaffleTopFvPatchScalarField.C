// Modified thermalBaffleTopFvPatchScalarField.C to allow TOP coupling (instead of default bottom)

#include "thermalBaffleTopFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyPolyPatch.H"
#include "mappedWallPolyPatch.H"

namespace Foam
{
namespace compressible
{

thermalBaffleTopFvPatchScalarField::thermalBaffleTopFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF),
    owner_(false),
    internal_(true),
    baffle_(nullptr),
    dict_(),
    extrudeMeshPtr_()
{}

thermalBaffleTopFvPatchScalarField::thermalBaffleTopFvPatchScalarField
(
    const thermalBaffleTopFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(ptf, p, iF, mapper),
    owner_(ptf.owner_),
    internal_(ptf.internal_),
    baffle_(nullptr),
    dict_(ptf.dict_),
    extrudeMeshPtr_()
{}

thermalBaffleTopFvPatchScalarField::thermalBaffleTopFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF, dict),
    owner_(false),
    internal_(true),
    baffle_(nullptr),
    dict_(dict),
    extrudeMeshPtr_()
{
    const fvMesh& thisMesh = patch().boundaryMesh().mesh();

    word regionName("none");
    dict_.readIfPresent("region", regionName);
    dict_.readIfPresent("internal", internal_);

    const word baffleName("3DBaffle" + regionName);

    if (!thisMesh.time().foundObject<fvMesh>(regionName) && regionName != "none")
    {
        if (!extrudeMeshPtr_)
        {
            createPatchMesh();
        }

        baffle_.reset(baffleType::New(thisMesh, dict_));
        owner_ = true;
        baffle_->rename(baffleName);
    }
}

thermalBaffleTopFvPatchScalarField::thermalBaffleTopFvPatchScalarField
(
    const thermalBaffleTopFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(ptf, iF),
    owner_(ptf.owner_),
    internal_(ptf.internal_),
    baffle_(nullptr),
    dict_(ptf.dict_),
    extrudeMeshPtr_()
{}

void thermalBaffleTopFvPatchScalarField::createPatchMesh()
{
    const fvMesh& thisMesh = patch().boundaryMesh().mesh();
    const word regionName(dict_.get<word>("region"));

    polyPatchList regionPatches(3);
    List<dictionary> dicts(regionPatches.size());
    List<word> patchNames(regionPatches.size());
    List<word> patchTypes(regionPatches.size());

    patchNames[bottomPatchID] = word("bottom");
    patchNames[sidePatchID] = word("side");
    patchNames[topPatchID] = word("top");

    patchTypes[topPatchID] = mappedWallPolyPatch::typeName; // TOP is coupled now

    if (internal_)
    {
        patchTypes[bottomPatchID] = mappedWallPolyPatch::typeName; // OPTIONAL
    }
    else
    {
        patchTypes[bottomPatchID] = polyPatch::typeName;
    }

    if (dict_.get<bool>("columnCells"))
    {
        patchTypes[sidePatchID] = emptyPolyPatch::typeName;
    }
    else
    {
        patchTypes[sidePatchID] = polyPatch::typeName;
    }

    const auto& mpp = refCast<const mappedPatchBase>(patch().patch(), dict_);
    const word coupleGroup(mpp.coupleGroup());

    wordList inGroups(1);
    inGroups[0] = coupleGroup;

    // NOW: TOP is coupled with fluid patch
    dicts[topPatchID].add("coupleGroup", coupleGroup);
    dicts[topPatchID].add("inGroups", inGroups);
    dicts[topPatchID].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);
    dicts[topPatchID].add("samplePatch", patch().name());
    dicts[topPatchID].add("sampleRegion", thisMesh.name());

    // Bottom: optionally coupled (if internal_) or free surface
    if (internal_)
    {
        const word coupleGroupSlave = coupleGroup + "_slave";
        inGroups[0] = coupleGroupSlave;
        dicts[bottomPatchID].add("coupleGroup", coupleGroupSlave);
        dicts[bottomPatchID].add("inGroups", inGroups);
        dicts[bottomPatchID].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);
    }

    forAll(regionPatches, patchi)
    {
        dictionary& patchDict = dicts[patchi];
        patchDict.set("nFaces", 0);
        patchDict.set("startFace", 0);

        regionPatches.set
        (
            patchi,
            polyPatch::New
            (
                patchTypes[patchi],
                patchNames[patchi],
                dicts[patchi],
                patchi,
                thisMesh.boundaryMesh()
            )
        );
    }

    extrudeMeshPtr_.reset
    (
        new extrudePatchMesh(thisMesh, patch(), dict_, regionName, regionPatches)
    );
}

void thermalBaffleTopFvPatchScalarField::updateCoeffs()
{
    if (this->updated()) return;
    if (owner_) baffle_->evolve();
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs();
}

void thermalBaffleTopFvPatchScalarField::write(Ostream& os) const
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::write(os);

    if (owner_)
    {
        os.writeEntry("extrudeModel", dict_.get<word>("extrudeModel"));
        os.writeEntry("nLayers", dict_.get<label>("nLayers"));
        os.writeEntry("expansionRatio", dict_.get<scalar>("expansionRatio"));
        os.writeEntry("columnCells", dict_.get<Switch>("columnCells"));

        const word extrudeModel(dict_.get<word>("extrudeModel") + "Coeffs");
        dict_.subDict(extrudeModel).writeEntry(extrudeModel, os);

        os.writeEntry("region", dict_.get<word>("region"));
        os.writeEntryIfDifferent<bool>("internal", true, internal_);
        os.writeEntry("active", dict_.get<Switch>("active"));

        dict_.subDict("thermoType").writeEntry("thermoType", os);
        dict_.subDict("mixture").writeEntry("mixture", os);
        dict_.subDict("radiation").writeEntry("radiation", os);
    }
}

makePatchTypeField(fvPatchScalarField, thermalBaffleTopFvPatchScalarField);

} // End namespace compressible
} // End namespace Foam

