/** \file Interpret.c

Contains the functions required to read the parameter file, and
functions that provide runtime information. The function var()
associates a string to a global variable. The function ReadVariables()
reads the content of a parameter file.  In addition, this file
contains a function that prints the command line usage to the standard
output, a function that provides verbose information about the setup
(if the -v switch is set on the command line), and functions that act
as a chronometer (if the -t switch is set on the command line).
*/

#include "mp.h"
#define MAXVARIABLES 500

extern int      begin_i;
extern boolean  OpenInner, OpenInnerDust;
static Param    VariableSet[MAXVARIABLES];
static int      VariableIndex = 0;
static int	FirstStep = YES;
static clock_t  First, Preceeding, Current, FirstUser, CurrentUser, PreceedingUser;
static long	Ticks;
boolean         FastTransport = YES, GuidingCenter = NO, BinaryCenter = NO, Indirect_Term = YES, Discard_GasIndirect_term = NO, RetrogradeBinary = NO;
boolean         IsDisk = YES, NonReflecting = NO, Corotating = NO, OuterSourceMass = NO, Evanescent = NO, MixedBC = NO, AccBoundary = NO;
boolean         Write_Density = YES, Write_Velocity = YES, Write_Energy = NO;
boolean         Write_Temperature = NO, Write_DivV = NO, Write_Jacobi = NO, Write_DustSystem = YES;
boolean         Write_TherDiff = NO, Write_RadDiff = NO, Write_TherCool = NO, Write_ViscHeat = NO, Write_DustDensity = NO;
boolean         Write_Potential = NO, Write_Test = NO, Write_gr = NO, Write_gtheta = NO;
boolean         Write_RadFBAcc = NO, Write_AziFBAcc = NO;
boolean         SelfGravity = NO, SGZeroMode = NO, ZMPlus = NO, AddNoise = NO;
boolean         EnergyEquation = NO, EntropyDiffusion = NO, RadiativeDiffusion = NO, ImplicitRadiativeDiffusion = NO, ThermalCooling = NO, StellarIrradiation = NO, ImposedStellarLuminosity = NO, ViscousHeating = YES, TempPresc = NO, BetaCooling = NO, SetConstantOpacity = NO;
boolean         CICPlanet = NO, ForcedCircular = NO, ForcedInnerCircular = NO, ComputeCPDMass = NO;
boolean         MHDLSA = NO, HighMCutoff = NO;
boolean         KNOpen = NO;
boolean         AdvecteLabel = NO;
boolean         SoftWriting = NO;
boolean         RetrogradePlanet = NO, AddMass = NO, ImposedDensity = NO;
boolean         ReadPlanetFileAtRestart = YES, DontApplySubKeplerian = NO;
boolean         NGPInterpolation = NO, CICInterpolation = NO, TSCInterpolation = YES;
boolean         DampToIni = NO, DampToAxi = YES, DampToViscous = NO;
boolean         CorotateWithOuterPlanet = NO;
boolean         DiscEvaporation = NO, ShortFrictionTimeApproximation = YES;
boolean         CustomizedIT = NO, AddFloors = YES, AddM1 = NO, AddM1Boosted = NO, AddM1toM10 = NO, TailOffGauss = NO, ExponentialCutoff = NO;
boolean         DustFeelDisk = YES, DustFeelSG = YES, DustFeelSGZeroMode = NO, DustFeelPlanets = YES, RestartWithNewDust = NO, DustFeelTurb = NO, TailOffIn = NO, TailOffAurelien = NO, TailOffStype = NO, TailOffGI = NO, TailOffSareh = NO, DustGrowth = NO;
boolean         ZZIntegrator = NO, NoTimestepConstraintByParticles = NO, PhotoEvaporation = NO, DecInner = NO, DustFeedback = NO, RemoveDustFromPlanetsHillRadius = YES;
boolean         DustFluid = NO, DustDiffusion = NO, Write_StokesNumber = NO;
boolean         BC1D_SS_ZeroVel = NO, BC1D_SS_NonZeroVel = NO, BC1D_ZeroDens = YES;
boolean         BM08 = NO;
boolean         CavityTorque = NO, CompareSGAndSummationTorques = NO;
boolean         ImposeZeroRadialVelocityBC = NO;

void var(name, ptr, type, necessary, deflt)
     char           *name;
     char           *ptr;
     int             type;
     int             necessary;
     char           *deflt;
{
  real            valuer;
  int             valuei;
  double	  temp;
  sscanf (deflt, "%lf", &temp);
  valuer = (real) (temp);
  valuei = (int) valuer;
  strcpy(VariableSet[VariableIndex].name, name);
  VariableSet[VariableIndex].variable = ptr;
  VariableSet[VariableIndex].type = type;
  VariableSet[VariableIndex].necessary = necessary;
  VariableSet[VariableIndex].read = NO;
  if (necessary == NO) {
    if (type == INT) {
      *((int *) ptr) = valuei;
    } else if (type == REAL) {
      *((real *) ptr) = valuer;
    } else if (type == STRING) {
      strcpy (ptr, deflt);
    }
  }
  VariableIndex++;
}

void ReadVariables(filename)
     char *filename;
{
  char            nm[300], s[350],stringval[290];
  char           *s1;
  double          temp;
  real            valuer;
  int             i, found, valuei, success, type;
  int            *ptri;
  real           *ptrr;
  FILE           *input;
  
  InitVariables();
  input = fopen(filename, "r");
  if (input == NULL) {
    mastererr ("Unable to read '%s'. Program stopped.\n",filename);
    prs_exit(1);
  }
  mastererr ("Reading parameters file '%s'.\n", filename);
  while (fgets(s, 349, input) != NULL) {
    success = sscanf(s, "%s ", nm);
    if ((nm[0] != '#') && (success == 1)) {  /* # begins a comment line */
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%lf", &temp);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%289s ", stringval);
      valuer = (real) temp;
      valuei = (int) temp;
      for (i = 0; i < strlen(nm); i++) {
	nm[i] = (char) toupper(nm[i]);
      }
      found = NO;
      for (i = 0; i < VariableIndex; i++) {
	if (strcmp(nm, VariableSet[i].name) == 0) {
	  if (VariableSet[i].read == YES) {
	    mastererr("Warning : %s defined more than once.\n", nm);
	  }
	  found = YES;
	  VariableSet[i].read = YES;
	  ptri = (int *) (VariableSet[i].variable);
	  ptrr = (real *) (VariableSet[i].variable);
	  if (VariableSet[i].type == INT) {
	    *ptri = valuei;
	  } else if (VariableSet[i].type == REAL) {
	    *ptrr = valuer;
	  } else if (VariableSet[i].type == STRING) {
	    strcpy (VariableSet[i].variable, stringval);
	  }
	}
      }
      if (found == NO) {
	mastererr("Warning : variable %s defined but non-existent in code.\n", nm);
      }
    }
  }
  
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ((VariableSet[i].read == NO) && (VariableSet[i].necessary == YES)) {
      if (found == NO) {
	mastererr("Fatal error : undefined mandatory variable(s):\n");
	found = YES;
      }
      mastererr("%s\n", VariableSet[i].name);
    }
    if (found == YES)
      prs_exit(1);
    
  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if (VariableSet[i].read == NO) {
      if (found == NO) {
	mastererr("Secondary variables omitted :\n");
	found = YES;
      }
      if ((type = VariableSet[i].type) == REAL)
	mastererr("%s ;\t Default Value : %.5g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
	mastererr("%s ;\t Default Value : %d\n", VariableSet[i].name, *((int *) VariableSet[i].variable));
      if (type == STRING)
	mastererr("%s ;\t Default Value : %s\n", VariableSet[i].name, VariableSet[i].variable);
    }
  }
  if ((*ADVLABEL == 'y') || (*ADVLABEL == 'Y')) AdvecteLabel = YES;
  if ((*OUTERSOURCEMASS == 'y') || (*OUTERSOURCEMASS == 'Y')) OuterSourceMass = YES;
  if ((*TRANSPORT == 's') || (*TRANSPORT == 'S')) FastTransport = NO;
  if ((*OPENINNERBOUNDARY == 'O') || (*OPENINNERBOUNDARY == 'o')) OpenInner = YES;
  if ((*OPENINNERBOUNDARY == 'N') || (*OPENINNERBOUNDARY == 'n')) NonReflecting = YES;
  if ((*OPENINNERBOUNDARY == 'E') || (*OPENINNERBOUNDARY == 'e')) Evanescent = YES;
  if ((*OPENINNERBOUNDARY == 'M') || (*OPENINNERBOUNDARY == 'm')) MixedBC = YES;
  if ((*OPENINNERBOUNDARY == 'A') || (*OPENINNERBOUNDARY == 'a')) AccBoundary = YES;
  if ((*OPENINNERBOUNDARY == 'K') || (*OPENINNERBOUNDARY == 'k')) {
    KNOpen = YES;
    OpenInner = YES;
  }
  if ((*OPENINNERBOUNDARYDUST == 'O') || (*OPENINNERBOUNDARYDUST == 'o')) OpenInnerDust = YES;
  // Boundary conditions 1D grid
  if ((*BOUNDARY1DGRID == 'Z') || (*BOUNDARY1DGRID == 'z')) {
    BC1D_SS_ZeroVel = YES;
    BC1D_ZeroDens = NO;
  }
  if ((*BOUNDARY1DGRID == 'D') || (*BOUNDARY1DGRID == 'd')) {
    BC1D_SS_NonZeroVel = YES;
    BC1D_ZeroDens = NO;
  }
  //
  if ((*TAILOFF == 'G') || (*TAILOFF == 'g')) {
    TailOffGauss = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*TAILOFF == 'H') || (*TAILOFF == 'h')) {
    TailOffSareh = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*TAILOFF == 'E') || (*TAILOFF == 'e')) {
    ExponentialCutoff = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*TAILOFF == 'A') || (*TAILOFF == 'a')) {
    TailOffAurelien = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*TAILOFF == 'B') || (*TAILOFF == 'b')) {
    TailOffStype = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*TAILOFF == 'S') || (*TAILOFF == 's')) {
    TailOffGI = YES;
    //CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*TAILOFF == 'I') || (*TAILOFF == 'i')) {
    TailOffIn = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*GRIDSPACING == 'L') || (*GRIDSPACING == 'l')) LogGrid = YES;
  if ((*DISK == 'N') || (*DISK == 'n')) IsDisk = NO;
  if ((*FRAME == 'C') || (*FRAME == 'c')) Corotating = YES;
  if ((*FRAME == 'G') || (*FRAME == 'g')) {
    Corotating = YES;
    GuidingCenter = YES;
  }
  if ((*FRAME == 'B') || (*FRAME == 'b')) {
    Corotating = YES;
    BinaryCenter = YES;
  }
  if ((*RETROGRADEBINARY == 'Y') || (*RETROGRADEBINARY == 'y')) RetrogradeBinary = YES;
  if ((*WRITEVELOCITY == 'N') || (*WRITEVELOCITY == 'n')) Write_Velocity = NO;
  if ((*WRITEDENSITY == 'N') || (*WRITEDENSITY == 'n')) Write_Density = NO;
  if ((*WRITEENERGY == 'Y') || (*WRITEENERGY == 'y')) Write_Energy = YES;
  if ((*WRITETEMPERATURE == 'Y') || (*WRITETEMPERATURE == 'y')) Write_Temperature = YES;
  if ((*WRITEDIVV == 'Y') || (*WRITEDIVV == 'y')) Write_DivV = YES;
  if ((*WRITEVISCHEAT == 'Y') || (*WRITEVISCHEAT == 'y')) Write_ViscHeat = YES;
  if ((*WRITETHERDIFF == 'Y') || (*WRITETHERDIFF == 'y')) Write_TherDiff = YES;
  if ((*WRITERADDIFF == 'Y') || (*WRITERADDIFF == 'y')) Write_RadDiff = YES;
  if ((*WRITETHERCOOL == 'Y') || (*WRITETHERCOOL == 'y')) Write_TherCool = YES;
  if ((*WRITEPOTENTIAL == 'Y') || (*WRITEPOTENTIAL == 'y')) Write_Potential = YES;
  if ((*WRITETEST == 'Y') || (*WRITETEST == 'y')) Write_Test = YES;
  if ((*WRITEGR == 'Y') || (*WRITEGR == 'y')) Write_gr = YES;
  if ((*WRITEGTHETA == 'Y') || (*WRITEGTHETA == 'y')) Write_gtheta = YES;
  if ((*WRITEJACOBI == 'Y') || (*WRITEJACOBI == 'y')) Write_Jacobi = YES;
  if ((*WRITEDUSTSYSTEM == 'N') || (*WRITEDUSTSYSTEM == 'n')) Write_DustSystem = NO;
  if ((*WRITEDUSTSTOKES == 'Y') || (*WRITEDUSTSTOKES == 'y')) Write_StokesNumber = YES;
  if ((*WRITERADFBACC == 'Y') || (*WRITERADFBACC == 'y')) Write_RadFBAcc = YES;
  if ((*WRITEAZIFBACC == 'Y') || (*WRITEAZIFBACC == 'y')) Write_AziFBAcc = YES;
  if ((*DUSTFLUID == 'Y') || (*DUSTFLUID == 'y')) DustFluid = YES;
  if (DustFluid)
    Write_DustDensity = YES;
  if ( (NBPART != 0) || (DustFluid) ) {
    Write_DustDensity = YES;
    if ((*WRITEDUSTDENSITY == 'N') || (*WRITEDUSTDENSITY == 'n')) Write_DustDensity = NO;
    /* Boolean for restarting a pre-evolved gas disc simulation with
       inclusion of dust particles */
    if ((*RESTARTWITHNEWDUST == 'y') || (*RESTARTWITHNEWDUST == 'Y')) RestartWithNewDust = YES;
  }
  if ((*DUSTDIFFUSION == 'Y') || (*DUSTDIFFUSION == 'y')) DustDiffusion = YES;
  if ((*SOFTWRITING == 'y') || (*SOFTWRITING == 'Y')) SoftWriting = YES;
  if ((*DUSTGROWTH == 'Y') || (*DUSTGROWTH == 'y')) DustGrowth = YES;
  if ((*DUSTFEELDISK == 'N') || (*DUSTFEELDISK == 'n')) DustFeelDisk = NO;
  if ((*REMOVEDUSTFROMPLANETSHILLRADIUS == 'N') || (*REMOVEDUSTFROMPLANETSHILLRADIUS == 'n')) RemoveDustFromPlanetsHillRadius = NO;
  if ((*DUSTFEELSG == 'N') || (*DUSTFEELSG == 'n')) DustFeelSG = NO;
  if ((*DUSTFEELSG == 'Z') || (*DUSTFEELSG == 'z')) {
    DustFeelSG = YES;
    DustFeelSGZeroMode = YES;
  }
  if ((*DUSTFEELPLANETS == 'N') || (*DUSTFEELPLANETS == 'n')) DustFeelPlanets = NO;
  if ((*DUSTFEELTURB == 'Y') || (*DUSTFEELTURB == 'y') || (*DUSTFEELTURB == 'C') || (*DUSTFEELTURB == 'c')) {
    DustFeelTurb = YES;
    srand48(time(NULL));
  }
  if ((*RETROGRADEPLANET == 'y') || (*RETROGRADEPLANET == 'Y')) RetrogradePlanet = YES;
  if ((*ADDMASS == 'y') || (*ADDMASS == 'Y')) AddMass = YES;
  if ((*DISCEVAPORATION == 'y') || (*DISCEVAPORATION == 'Y')) {
    DiscEvaporation = YES;
    masterprint("Disc evaporation included, with characteristic timescale = %lg\n", TEVAP);
  }
  if ((*CUSTIT == 'y') || (*CUSTIT == 'Y')) {
    CustomizedIT = YES;
  }
  if ((*ADDFLOORS == 'n') || (*ADDFLOORS == 'N')) {
    AddFloors = NO;
  }
  if ((*ADDM1 == 'y') || (*ADDM1 == 'Y')) {
    AddM1 = YES;
  }
  if ((*ADDM1 == 'b') || (*ADDM1 == 'B')) {
    AddM1Boosted = YES;
  }
  if ((*ADDM1TOM10 == 'y') || (*ADDM1TOM10 == 'Y')) {
    AddM1toM10 = YES;
  }
  if ((*ZZINTEGRATOR == 'y') || (*ZZINTEGRATOR == 'Y')) {
    ZZIntegrator = YES;
  }
  if ((*PHOTOEVAPORATION == 'Y') || (*PHOTOEVAPORATION == 'y')) {
    PhotoEvaporation = YES;
  }
  if ( (PhotoEvaporation == YES) && ( LX == 0.0)) {
    mastererr("You should set the star X-ray luminosity in erg/s for photoevaporation in the .par parameter file.\n");
    prs_exit(1);
  }
  if ((*DECINNER == 'y') || (*DECINNER == 'Y')) {
    DecInner = YES;
  }
  if ((*DUSTFEEDBACK == 'y') || (*DUSTFEEDBACK == 'Y')) {
    DustFeedback = YES;
  }
  if ((*INTERPOLATION == 'N') || (*INTERPOLATION == 'n')) {
    NGPInterpolation = YES;
    CICInterpolation = NO;
    TSCInterpolation = NO;
  }
  if ((*INTERPOLATION == 'C') || (*INTERPOLATION == 'c')) {
    NGPInterpolation = NO;
    CICInterpolation = YES;
    TSCInterpolation = NO;
  }
  if ((*INTERPOLATION == 'T') || (*INTERPOLATION == 't')) {
    NGPInterpolation = NO;
    CICInterpolation = NO;
    TSCInterpolation = YES;
  }
  if ((*INTERPOLATION == 'T') || (*INTERPOLATION == 't')) {
    NGPInterpolation = NO;
    CICInterpolation = NO;
    TSCInterpolation = YES;
  }
  if ((*NODTCONSTRAINTBYPCS == 'y') || (*NODTCONSTRAINTBYPCS == 'Y')) {
    NoTimestepConstraintByParticles = YES;
  }
  if ((*SFTAPPROX == 'n') || (*SFTAPPROX == 'N')) {
    ShortFrictionTimeApproximation = NO;
  }
  if ((*DAMPTOINI == 'y') || (*DAMPTOINI == 'Y')) {
    DampToIni = YES;
    DampToAxi = NO;
    DampToViscous = NO;
  }
  if ((*DAMPTOAXI == 'n') || (*DAMPTOAXI == 'N')) {
    DampToAxi = NO;
    DampToIni = YES;
    DampToViscous = NO;
  }
  if ((*DAMPTOVISCOUS == 'y') || (*DAMPTOVISCOUS == 'Y')) {
    DampToViscous = YES;
    DampToIni = NO;
    DampToAxi = NO;
  }
  if ((*COROTATEWITHOUTERPLANET == 'y') || (*COROTATEWITHOUTERPLANET == 'Y')) CorotateWithOuterPlanet = YES;
  if ((*INDIRECTTERM == 'N') || (*INDIRECTTERM == 'n')) Indirect_Term = NO;
  if ((*DISCARDGASINDIRECTTERM == 'Y') || (*DISCARDGASINDIRECTTERM == 'y')) Discard_GasIndirect_term = YES;
  if ((*SELFGRAVITY == 'Y') || (*SELFGRAVITY == 'y')) SelfGravity = YES;
  if ((*SELFGRAVITY == 'Z') || (*SELFGRAVITY == 'z')) {
    SelfGravity = YES;
    SGZeroMode = YES;
  }
  if ((*ZMPLUS == 'Y') || (*ZMPLUS == 'y')) ZMPlus = YES;
  if ( (ZMPlus) && (!SGZeroMode) ) {
    masterprint ("This is not very meaningfull to involve the anisotropic pressure model (ZMPlus=Yes) without taking into account the axisymmetric component of the disk self-gravity. I decided to put ZMPlus = No. Please check again!");
    ZMPlus = NO;
  }
  if ((*ADDNOISE == 'Y') || (*ADDNOISE == 'y')) {
    AddNoise = YES;
    srand48(time(NULL));
  }
  if ((*ENERGYEQUATION == 'Y') || (*ENERGYEQUATION == 'y')) {
    EnergyEquation = YES;
    Write_Temperature = YES;
  }
  if ((*ENTROPYDIFFUSION == 'Y') || (*ENTROPYDIFFUSION == 'y')) EntropyDiffusion = YES;
  if ((*RADIATIVEDIFFUSION == 'E') || (*RADIATIVEDIFFUSION == 'e')) {
    RadiativeDiffusion = YES;
  }
  if ((*RADIATIVEDIFFUSION == 'I') || (*RADIATIVEDIFFUSION == 'i')) {
    RadiativeDiffusion = YES;
    ImplicitRadiativeDiffusion = YES;
  }
  if ((*THERMALCOOLING == 'Y') || (*THERMALCOOLING == 'y')) ThermalCooling = YES;
  if ((*STELLARIRRADIATION == 'Y') || (*STELLARIRRADIATION == 'y')) StellarIrradiation = YES;
  if ((*STELLARIRRADIATION == 'L') || (*STELLARIRRADIATION == 'l')) {
    StellarIrradiation = YES;
    ImposedStellarLuminosity = YES;
  }
  if ((*SETCONSTANTOPACITY == 'Y') || (*SETCONSTANTOPACITY == 'y')) SetConstantOpacity = YES;
  if ((*VISCOUSHEATING == 'N') || (*VISCOUSHEATING == 'n')) ViscousHeating = NO;
  if ((*TEMPPRESC == 'Y') || (*TEMPPRESC == 'y')) TempPresc = YES;
  if ((*BETACOOLING == 'Y') || (*BETACOOLING == 'y')) BetaCooling = YES;
  if ( (EnergyEquation) && (ADIABATICINDEX == 1) ) {
    masterprint ("You cannot have EnergyEquation = YES and AdiabatcIndex = 1. I decided to put EnergyEquation = No, to simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
    EnergyEquation = NO;
  }
  if ((*WRITEENERGY == 'N') || (*WRITEENERGY == 'n')) Write_Energy = NO;
  if ((*MHD == 'L') || (*MHD == 'l')) {
    MHDLSA = YES;
    srand48(time(NULL));
  }
  if ((*HIGHMCUTOFF == 'Y') || (*HIGHMCUTOFF == 'y')) {
    HighMCutoff = YES;
  }
  if ((*BM08TRICK == 'Y') || (*BM08TRICK == 'y')) {
    BM08 = YES;
  }
  if ((*CAVITYTORQUE == 'Y') || (*CAVITYTORQUE == 'y')) {
    /* June 2022 (Guillaume Robert's internship) */
    CavityTorque = YES;
    CentrifugalBalance = YES;
    DontApplySubKeplerian = YES;
  }
  if ((*ZEROVRADBC == 'Y') || (*ZEROVRADBC == 'y')) {
    ImposeZeroRadialVelocityBC = YES;
  }
  if ((*COMPARESGANDSUMMATIONTORQUES == 'Y') || (*COMPARESGANDSUMMATIONTORQUES == 'y')) {
    CompareSGAndSummationTorques = YES;
  }
  if ((*IMPOSEDDENSITY == 'y') || (*IMPOSEDDENSITY == 'Y')) ImposedDensity = YES;
  if ((*EXCLUDEHILL == 'Y') || (*EXCLUDEHILL == 'y')) ExcludeHill = YES;
  if ((*CICPLANET == 'Y') || (*CICPLANET == 'y')) CICPlanet = YES;
  if ((*FORCEDCIRCULAR == 'Y') || (*FORCEDCIRCULAR == 'y')) ForcedCircular = YES;
  if ((*FORCEDINNERCIRCULAR == 'Y') || (*FORCEDINNERCIRCULAR == 'y')) ForcedInnerCircular = YES;
  if ((*COMPUTECPDMASS == 'Y') || (*COMPUTECPDMASS == 'y')) ComputeCPDMass = YES;
  if ((*READPLANETFILEATRESTART == 'N') || (*READPLANETFILEATRESTART == 'n')) ReadPlanetFileAtRestart = NO;
  if ((*DONTAPPLYSUBKEPLERIAN == 'Y') || (*DONTAPPLYSUBKEPLERIAN == 'y')) {
    masterprint ("================================================\n");
    masterprint ("I will not apply subkeplerian boundary on vtheta\n");
    masterprint ("================================================\n");
    DontApplySubKeplerian = YES;
  }
  if (Evanescent) {
    masterprint ("Evanescent wave-killing zones boundary condition is applied,\n");
    masterprint ("I will therefore not apply subKeplerian boundary condition on vtheta.\n");
    DontApplySubKeplerian = YES;
  }
  if ( (EXCLUDEHILLFACTOR < 0.0) || (EXCLUDEHILLFACTOR > 1.0) ) {
    mastererr ("EXCLUDEHILLFACTOR must range between 0 and 1.\n");
    prs_exit (1);
  }
  if ((ALPHAVISCOSITY != 0.0) && (VISCOSITY != 0.0)) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("VISCOSITY and ALPHAVISCOSITY.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  }
  if (ALPHAVISCOSITY != 0.0) {
    ViscosityAlpha = YES;
    masterprint ("Viscosity is of alpha type\n");
  }
  if (SelfGravity && (SGTHICKNESSSMOOTHING == 0.0)) {
    mastererr ("You cannot have a vanishing smoothing length for \n");
    mastererr ("the self-gravitating kernel. Please modify and run again.\n");
    prs_exit (1);
  }
  if (ROCHESMOOTHING != 0.0) {
    RocheSmoothing = YES;
    masterprint ("Planet potential smoothing scales with their Hill sphere.\n");
  }
  if (OverridesOutputdir == YES) {
    sprintf (OUTPUTDIR, "%s", NewOutputdir);
  }
  if (Evanescent && ( (WKZRMIN == 0.0) || (WKZRMAX == 0.0) )) {
    mastererr ("Evanescent 'wave-killing zones' assumed as a boundary condition, \n");
    mastererr ("but you did not specify the radii of the border zones. Please run again.\n");
    mastererr ("by setting the values for WKZRMIN and WKZRMAX");
    prs_exit (1);
  }
  if (MixedBC) {
    WKZRMIN = 0.0;
    if (WKZRMAX == 0.0) {
      mastererr ("Evanescent 'wave-killing zones' assumed as the outer boundary condition, \n");
      mastererr ("but you did not specify the radii of the border zones. Please run again.\n");
      mastererr ("by setting the value of WKZRMAX");
      prs_exit (1);
    }
  }
  /* Add a trailing slash to OUTPUTDIR if needed */
  if (*(OUTPUTDIR+strlen(OUTPUTDIR)-1) != '/')
    strcat (OUTPUTDIR, "/");
}

void PrintUsage (execname)
     char *execname;
{
  mastererr("Usage : %s [-abcdeimnptvz] [-(0-9)] [-s number] [-f scaling] parameters file\n", execname);
  mastererr("\n-a : Monitor mass and angular momentum at each timestep\n");
  mastererr("-b : Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n");
  mastererr("-c : Sloppy CFL condition (checked at each DT, not at each timestep)\n");
  mastererr("-d : Print some debugging information on 'stdout' at each timestep\n");
  mastererr("-e : Activate EU test problem torque file output\n");
  mastererr("-f : Scale density array by 'scaling'. Useful to increase/decrease\n");
  mastererr("     disk surface density after a restart, for instance.            \n");
  mastererr("-i : tabulate Sigma profile as given by restart files\n");
  mastererr("-m : Merge output files from different CPUs\n");
  mastererr("-n : Disable simulation. The program just reads parameters file\n");
  mastererr("-o : Overrides output directory of input file.\n");
  mastererr("-p : Give profiling information at each time step\n");
  mastererr("-s : Restart simulation, taking #'number' files as initial conditions\n");
  mastererr("-t : Monitor CPU time usage at each time step\n");
  mastererr("-v : Verbose mode. Tells everything about parameters file\n");
  mastererr("-z : fake sequential built when evaluating sums on HD meshes\n");
  mastererr("-(0-9) : only write initial (or restart) HD meshes,\n");
  mastererr("     proceed to the next nth output and exit\n");
  mastererr("     This option must stand alone on one switch (-va -4 is legal, -v4a is not)\n");
  prs_exit (1);
}

real TellNbOrbits (time)
     real time;
{
  return time/(2.0*M_PI)*sqrt(G*1.0/1.0/1.0/1.0);
}

real TellNbOutputs (time)
     real time;
{
  return (time/DT/NINTERM);
}

void TellEverything () {
  real temp, nbfileoutput;
  if (!CPU_Master) return;
  printf ("\nDisc properties:\n");
  printf ("----------------\n");
  printf ("Inner Radius          : %g\n", RMIN);
  printf ("Outer Radius          : %g\n", RMAX);
  printf ("Aspect Ratio          : %g\n", ASPECTRATIO);
  printf ("VKep at inner edge    : %.3g\n", sqrt(G*1.0*(1.-0.0)/RMIN));
  printf ("VKep at outer edge    : %.3g\n", sqrt(G*1.0/RMAX));
  temp=(PMAX-PMIN)*SIGMA0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - pow(RMIN,2.0-SIGMASLOPE));	/* correct this and what follows... */
  printf ("Initial Disk Mass             : %g\n", temp);
  temp=(PMAX-PMIN)*SIGMA0/(2.0-SIGMASLOPE)*(1.0 - pow(RMIN,2.0-SIGMASLOPE));
  printf ("Initial Mass inner to r=1.0  : %g \n", temp);
  temp=(PMAX-PMIN)*SIGMA0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - 1.0);
  printf ("Initial Mass outer to r=1.0  : %g \n", temp);
  printf ("Travelling time for acoustic density waves :\n");
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(RMIN,1.5));
  printf (" * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(1.0,1.5));
  printf (" * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(1.0,1.5)-pow(RMIN,1.5));
  printf (" * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = (PMAX-PMIN)*sqrt(RMIN*RMIN*RMIN/G/1.0);
  printf ("Orbital time at Rmin  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  temp = (PMAX-PMIN)*sqrt(RMAX*RMAX*RMAX/G/1.0);
  printf ("Orbital time at Rmax  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  printf ("Sound speed :\n");
  printf (" * At unit radius     : %.3g\n", ASPECTRATIO*sqrt(G*1.0));
  printf (" * At outer edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMAX));
  printf (" * At inner edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMIN));
  printf ("\nGrid properties:\n");
  printf ("----------------\n");
  printf ("Number of rings       : %d\n", NRAD);
  printf ("Number of sectors     : %d\n", NSEC);
  printf ("Total cells           : %d\n", NRAD*NSEC);
  printf ("\nOutputs properties:\n");
  printf ("-------------------\n");
  printf ("Time increment between outputs : %.3f = %.3f orbits\n", NINTERM*DT, TellNbOrbits(NINTERM*DT));
  printf ("At each output #i, the following files are written:\n");
  printf ("gasdens[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("gasvrad[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("gasvtheta[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  if (EnergyEquation == YES)
    printf ("Temperature[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  if (AdvecteLabel == YES)
    printf ("gaslabel[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("There will be in total %d outputs\n", NTOT/NINTERM);
  printf ("(which correspond to an elapsed time = %.3f or to %.2f orbits)\n", NTOT*DT, TellNbOrbits(NTOT*DT));
  nbfileoutput = 3.0;
  if (EnergyEquation == YES)
    nbfileoutput += 1.0;
  if (AdvecteLabel == YES)
    nbfileoutput += 1.0;
  temp =nbfileoutput*GLOBALNRAD*NSEC*sizeof(real);
  temp *= (real)NTOT/(real)NINTERM;
  temp /= 1024.0*1024.0;
  printf ("So the code will produce ~%.2f Mbytes of data\n", temp);
  printf ("Check (eg by issuing a 'df' command) that you have enough disk space,\n");
  printf ("otherwise you will get a system full and the code will stop.\n");
  fflush (stdout);
}

void GiveTimeInfo (number)
     int number;
{
  struct tms buffer;
  real total, last, mean, totalu;
  Current = times (&buffer);
  CurrentUser = buffer.tms_utime;
  if (FirstStep == YES) {
    First = Current;
    FirstUser = CurrentUser;
    fprintf (stderr, "Time counters initialized\n");
    FirstStep = NO;
    Ticks = sysconf (_SC_CLK_TCK);
  }
  else {
    total = (real)(Current - First)/Ticks;
    totalu= (real)(CurrentUser-FirstUser)/Ticks;
    last  = (real)(CurrentUser - PreceedingUser)/Ticks;
    number -= begin_i/NINTERM;
    mean  = totalu / number;
    fprintf (stderr, "Total Real Time elapsed    : %.3f s\n", total);
    fprintf (stderr, "Total CPU Time of process  : %.3f s (%.1f %%)\n", totalu, 100.*totalu/total);
    fprintf (stderr, "CPU Time since last time step : %.3f s\n", last);
    fprintf (stderr, "Mean CPU Time between time steps : %.3f s\n", mean);
    fprintf (stderr, "CPU Load on last time step : %.1f %% \n", (real)(CurrentUser-PreceedingUser)/(real)(Current-Preceeding)*100.);
    
  }	
  PreceedingUser = CurrentUser;
  Preceeding = Current;
}

void InitSpecificTime (profiling, process_name, title)
     boolean profiling;
     TimeProcess *process_name;
     char *title;
{
  struct tms buffer;
  if (profiling == NO) return;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  process_name->clicks = buffer.tms_utime;
  strcpy (process_name->name, title);
}

void GiveSpecificTime (profiling, process_name)
     boolean profiling;
     TimeProcess process_name;
{
  struct tms buffer;
  long ticks;
  real t;
  if (profiling == NO) return;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  ticks = buffer.tms_utime - process_name.clicks;
  t = (real)ticks / (real)Ticks;
  fprintf (stderr, "Time spent in %s : %.3f s\n", process_name.name, t);
}

