///////////////////////////////////////////////////////////////
////////////////////// TIMESTEP CONTROL ///////////////////////
///////////////////////////////////////////////////////////////

word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0;
else if (CFL == "Todd") CFLint = 1;
else if (CFL == "Courant") CFLint = 2; 
else
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "CFL condition unrecongnized : Coats, Todd and Courant available" 
            << exit(FatalError);
}

//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << "Reading pressure field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << nl << "Reading saturation field Sb" << endl;
volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));

scalar minDeltaT(runTime.controlDict().lookupOrDefault<scalar>("minDeltaT",0.));

//- Optional steady-state detection: if > 0, the run stops as soon as the
//  max change of Sb between two steps drops below this value (disabled
//  by default, opt-in via "convergence" in controlDict).
scalar convergence(runTime.controlDict().lookupOrDefault<scalar>("convergence",0.));
volScalarField SbOld
(
    IOobject
    (
        "SbOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Sb
);

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- list that receives event files of event-based boundary conditions
List<patchEventFile*> patchEventList;
eventInfiltration::setEventFileRegistry(&patchEventList, Sb.name());

dimensionedScalar Sbmin(transportProperties.lookupOrDefault("Sbmin",dimensionedScalar("Sbmin",dimless,0)));

//- Optional hard constraint: reset Sb to a fixed value in a given cellZone
//  every timestep. A fixedValue boundary condition alone only pins the
//  boundary face, not the adjacent cell centres, so it cannot be used to
//  represent e.g. a permanently saturated reservoir region; this clamp
//  operates directly on the cell values instead.
Switch fixedSbActive(false);
label fixedSbZoneID(-1);
scalar fixedSbValue(0);
if (transportProperties.found("fixedSb"))
{
    const dictionary& fixedSbDict = transportProperties.subDict("fixedSb");
    fixedSbActive = fixedSbDict.lookupOrDefault<Switch>("active", false);
    if (fixedSbActive)
    {
        word fixedSbZoneName(fixedSbDict.lookup("zoneName"));
        fixedSbZoneID = mesh.cellZones().findZoneID(fixedSbZoneName);
        if (fixedSbZoneID == -1)
        {
            FatalErrorIn("createFields.H")
                << "cellZone " << fixedSbZoneName
                << " not found for fixedSb constraint" << exit(FatalError);
        }
        fixedSbValue = readScalar(fixedSbDict.lookup("value"));
        Info << nl << "Clamping Sb = " << fixedSbValue << " in cellZone "
            << fixedSbZoneName << " every timestep" << endl;

        forAll(mesh.cellZones()[fixedSbZoneID], i)
        {
            Sb[mesh.cellZones()[fixedSbZoneID][i]] = fixedSbValue;
        }
        Sb.correctBoundaryConditions();
    }
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(mesh,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(mesh,transportProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();

dimensionedScalar Mmu(mub/mua);//viscosity ratio

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Relative permeability / Capillary pressure models 
autoPtr<twophasePorousMediumModel> porousModel = twophasePorousMediumModel::New("Sb", mesh, transportProperties, phaseb);
autoPtr<relativePermeabilityModel>& krModel = porousModel->krModel();
autoPtr<capillarityModel>& pcModel = porousModel->pcModel();
scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));

// Porosity
porousModel->check_eps();
const volScalarField& eps = porousModel->eps();

// Intrinsic permeability
porousModel->check_K();
const volScalarField& K = porousModel->K();

// permeability interpolation
surfaceScalarField Kf(fvc::interpolate(K,"K"));

/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info << nl << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ua + Ub
);

#include "createPhi.H"
surfaceScalarField phiP = phi;

///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless/dimTime,0)
);

/////////////////////////////////////////////////////////////////////////////
///////////////// READABLE LIQUID CONTENT DIAGNOSTIC FIELDS //////////////////
/////////////////////////////////////////////////////////////////////////////

//- Volumetric water content theta = eps*Sb [m3 liquid / m3 total], the
//  standard porous-media "moisture content" - easier to read than the raw
//  saturation Sb since it accounts for how much pore space actually exists.
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    eps*Sb
);

//- Mass of liquid phase b per cell [kg]
volScalarField mLiquid
(
    IOobject
    (
        "mLiquid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mLiquid", dimMass, 0)
);
mLiquid.primitiveFieldRef() = rhob.value()*theta.primitiveField()*mesh.V();

//- Running log of total liquid volume/mass in the whole domain, appended
//  at each write time (units chosen for a lab-scale sample: mL and g).
autoPtr<OFstream> liquidBalanceFile;
if (Pstream::master())
{
    liquidBalanceFile.reset(new OFstream("liquidBalance.csv"));
    liquidBalanceFile() << "time,totalLiquidVolume_mL,totalLiquidMass_g" << endl;
}
