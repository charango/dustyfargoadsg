/** \file SourceEuler.c 

    Contains routines used by the hydrodynamical loop. More specifically,
    it contains the main loop itself and all the source term substeps
    (with the exception of the evaluation of the viscous force). The
    transport substep is treated elsewhere. */

#include "mp.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

#define BINARYSECURITY 100.0

#define MAXITS 100000  /* maximum number of iteration in SOR solver for
			implicite radiative diffusion */

static PolarGrid *TemperInt;
static PolarGrid *VradNew, *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid *DTemperInt;
static PolarGrid *DVradNew, *DVradInt;
static PolarGrid *DVthetaNew, *DVthetaInt;
static PolarGrid *EnergyNew, *EnergyInt, *TempInt;

extern boolean Corotating;
extern boolean EnergyEquation, EntropyDiffusion, RadiativeDiffusion, ImplicitRadiativeDiffusion, ThermalCooling, ViscousHeating;
extern boolean SelfGravity, ZMPlus;
extern boolean DustDiffusion;
extern boolean FastTransport, IsDisk, AddFloors, DustFluid, DustFeedback, DustFeelDisk, DampToViscous, ShortFrictionTimeApproximation;
extern boolean AddMass, DiscEvaporation, PhotoEvaporation, CavityTorque;

Pair DiskOnPrimaryAcceleration;

int FirstGasStepFLAG=1;
static int AlreadyCrashed = 0, GasTimeStepsCFL;
static int niterSOR;
static int niterbuf = MAXITS;
static int jchess1st, jchess2nd;

static real timeCRASH;  
static real omegaSOR = 1.97;   // initial guess value for omega relaxation parameter in SOR solver
static real epsilonSOR = 1e-4; // relatice accuracy required for SOR solver


boolean DetectCrash (array)
     PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0) 
	bool = YES;
    }
  }
  return bool;
}
 
void FillPolar1DArrays ()
{
  extern boolean Restart;
  FILE *input, *output;
  int i, ii, j;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256], OutputName2[256], command[1024];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  if (!Restart) 
    sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  else
    sprintf (OutputName, "%s%s", OUTPUTDIR, "newused_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = 0.5*(PMAX-PMIN)*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  /* NEW (Feb. 2015): we output min and max radii handled by each
     CPU */
  sprintf (OutputName2, "%sminmaxradii.dat.%05d",OUTPUTDIR,CPU_Rank);
  output = fopen (OutputName2, "w");
  if (output == NULL) {
    mastererr ("Can't write %s.\nProgram stopped.\n", OutputName2);
    prs_exit (1);
  }
  fprintf (output, "%d\t%.18g\t%.18g\n", CPU_Rank, Rinf[Zero_or_active], Rsup[Max_or_active-1]);
  fclose (output);
  MPI_Barrier (MPI_COMM_WORLD);
  if (CPU_Master) {
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; sync; cat minmaxradii.dat.%05d >> minmaxradii.dat.%05d", \
	       OUTPUTDIR, i, CPU_Rank);
      system (command);
    }
    sprintf (command, "cd %s; mv minmaxradii.dat.%05d minmaxradii.dat; rm -f minmaxradii.dat.0*", \
    	     OUTPUTDIR, CPU_Rank);
    system (command);
  }
  /* -------------- */
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
  /* ---------------------------------- */
  /* New (Apr 12 2010): definition of a global array with azimuths of
     cells centres */
  /* ---------------------------------- */
  for (j = 0; j < NSEC; j++) {
    Azimuth[j] = PMIN + (PMAX-PMIN)*(real)j/(real)NSEC;
    CosAzimuth[j] = cos(Azimuth[j]);
    SinAzimuth[j] = sin(Azimuth[j]);
    /* case where azimuthal extent smaller than 2pi */
    //if ( fabs(PMAX-PMIN-2.*M_PI) > 0.2 ) 
    // Azimuth[j] = PMIN + (PMAX-PMIN)*(real)j/(real)(NSEC-1);
  }
  /* Global arrays with azimuths of cells interfaces */
  for (j = 0; j < NSEC; j++) {
    AziInf[j] = Azimuth[j]-0.5*(PMAX-PMIN)/(real)NSEC;
    AziSup[j] = Azimuth[j]+0.5*(PMAX-PMIN)/(real)NSEC;
  }
  
  if (CPU_Master) {
    sprintf (OutputName, "%s%s", OUTPUTDIR, "used_azi.dat");
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (j = 0; j < NSEC; j++) {
      fprintf (output, "%.18g\t%.18g\t%.18g\n", Azimuth[j],AziInf[j],AziSup[j]);
    }
    fclose (output);
  }
}

void InitEuler (Vr, Vt, Rho, Energy, DVr, DVt, DRho, sys)
     PolarGrid *Vr, *Vt, *Rho, *Energy;
     PolarGrid *DVr, *DVt, *DRho;
     PlanetarySystem *sys;
{
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  RadiativeKCoeff = CreatePolarGrid(NRAD, NSEC, "RadiativeKCoeff");
  RadiativeChiCoeff = CreatePolarGrid(NRAD, NSEC, "RadiativeChiCoeff");
  TempInt      = CreatePolarGrid(NRAD, NSEC, "TempInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  IndPotential = CreatePolarGrid(NRAD, NSEC, "IndirectPotential");
  RadGradP     = CreatePolarGrid(NRAD, NSEC, "RadGradP");
  AziGradP     = CreatePolarGrid(NRAD, NSEC, "AziGradP");
  RadIndAcc    = CreatePolarGrid(NRAD, NSEC, "RadIndAcc");
  AziIndAcc    = CreatePolarGrid(NRAD, NSEC, "AziIndAcc");
  RadSGAcc     = CreatePolarGrid(NRAD, NSEC, "RadSGAcc");
  AziSGAcc     = CreatePolarGrid(NRAD, NSEC, "AziSGAcc");
  RadFBAcc     = CreatePolarGrid(NRAD, NSEC, "RadFBAcc");
  AziFBAcc     = CreatePolarGrid(NRAD, NSEC, "AziFBAcc");
  FBedot       = CreatePolarGrid(NRAD, NSEC, "FBedot");
  TurbPotential= CreatePolarGrid(NRAD, NSEC, "TurbPotential");
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
  ViscHeat     = CreatePolarGrid(NRAD, NSEC, "ViscousHeating");
  EntropyDiff  = CreatePolarGrid(NRAD, NSEC, "EntropyDiff");
  RadiativeDiff= CreatePolarGrid(NRAD, NSEC, "RadiativeDiff");
  ThermCool    = CreatePolarGrid(NRAD, NSEC, "ThermalCooling");
  Opacity      = CreatePolarGrid(NRAD, NSEC, "Opacity");
  Test         = CreatePolarGrid(NRAD, NSEC, "Test");
  /* Case where dust is modelled as a low-pressure fluid */
  DRhoStar     = CreatePolarGrid(NRAD, NSEC, "DRhoStar"); /*initialize dust density star/int and vradint vradnew and vthetanew vthetaint */
  DRhoInt      = CreatePolarGrid(NRAD, NSEC, "DRhoInt");
  DVradNew     = CreatePolarGrid(NRAD, NSEC, "DVradNew");
  DVradInt     = CreatePolarGrid(NRAD, NSEC, "DVradInt");
  DVthetaNew   = CreatePolarGrid(NRAD, NSEC, "DVthetaNew");
  DVthetaInt   = CreatePolarGrid(NRAD, NSEC, "DVthetaInt");
  DTemperInt   = CreatePolarGrid(NRAD, NSEC, "DTemperInt");
  DPressure    = CreatePolarGrid(NRAD, NSEC, "DPressure");
  DSoundSpeed  = CreatePolarGrid(NRAD, NSEC, "DSoundSpeed");
  Stokes       = CreatePolarGrid(NRAD, NSEC, "Stokes"); /* dust stokes number */
  Diag1        = CreatePolarGrid(NRAD, NSEC, "Diag1"); /* three diagnostic */
  Diag2        = CreatePolarGrid(NRAD, NSEC, "Diag2");
  Diag3        = CreatePolarGrid(NRAD, NSEC, "Diag3");
  if (DustDiffusion == YES) {
    Fdiffrp         = CreatePolarGrid(NRAD, NSEC, "Fdiffrp");
    Fdifftp         = CreatePolarGrid(NRAD, NSEC, "Fdifftp");
  }
  if (RadiativeDiffusion && ImplicitRadiativeDiffusion) {
    a_SORarray      = CreatePolarGrid(NRAD, NSEC, "a_SORarray");
    b_SORarray      = CreatePolarGrid(NRAD, NSEC, "b_SORarray");
    c_SORarray      = CreatePolarGrid(NRAD, NSEC, "c_SORarray");
    d_SORarray      = CreatePolarGrid(NRAD, NSEC, "d_SORarray");
    e_SORarray      = CreatePolarGrid(NRAD, NSEC, "e_SORarray"); 
    f_SORarray      = CreatePolarGrid(NRAD, NSEC, "f_SORarray");
    Global_tempint         = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_tempint");
    Global_a_SORarray      = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_a_SORarray");
    Global_b_SORarray      = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_b_SORarray");
    Global_c_SORarray      = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_c_SORarray");
    Global_d_SORarray      = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_d_SORarray");
    Global_e_SORarray      = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_e_SORarray"); 
    Global_f_SORarray      = CreatePolarGrid(GLOBALNRAD, NSEC, "Global_f_SORarray");
  }
  InitComputeAccel ();
  /* Rho and Energy are already initialized: cf main.c */
  if (DustFluid) {
    ComputeDustSoundSpeedAndPressure (DRho);
    InitDustStokesNumber (Rho);
    InitDustVelocities (DVr, DVt, Rho, DRho);
  }
  /* We need to compute the sound speed at least once for locally
     isothermal runs */
  ComputeSoundSpeed (Rho, Energy, sys);
  ComputePressureField (Rho, Energy);
  ComputeTemperatureField (Rho, Energy);
  ComputeOpacities (Rho, Energy);
  InitGasVelocities (Vr, Vt, Rho, DRho);
  /* To output heating source terms at t=0... */
  ComputeViscousTerms (Vr, Vt, Rho, DVr, DVt, DRho);
  if (EnergyEquation && ViscousHeating)
    ComputeViscousHeating (Rho);
  if (EnergyEquation && EntropyDiffusion)
    ComputeEntropyDiffusion (Rho, Energy);
  if (EnergyEquation && RadiativeDiffusion && !ImplicitRadiativeDiffusion)
      ComputeRadiativeDiffusion (Rho, Energy);
  // CB (Dec 2017): function ComputeThermalCooling is no longer 
  // used (see substep 3, an implicit solver is used there)
  /*
  if (EnergyEquation && ThermalCooling)
    ComputeThermalCooling (Rho, Energy);
  */
}

real min2 (a,b)
     real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
     real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
     PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}

void AlgoGas (force, Rho, Vrad, Vtheta, Energy, Label, DRho, dustpcdens, DVrad, DVtheta, sys, dsys)
     Force *force;
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
     PolarGrid *DRho, *DVrad, *DVtheta, *dustpcdens;
     PlanetarySystem *sys;
     DustSystem *dsys;
{
  real dthydro, dtnbody, dt, buf, dtemp=0.0;
  real xk, xj, yk, yj, mk, mj, dist;
  real OmegaNew, domega;
  int gastimestepcfl, i, k, j, NbPlanets;
  boolean Crashed=NO;
  extern boolean ComputeCPDMass, NoTimestepConstraintByParticles, Evanescent, Write_DustDensity;
  real DustTimeStep;
  DustTimeStep = DUSTTIMESTEP;
  FirstGasStepFLAG=1;
  gastimestepcfl = 1;
  NbPlanets = sys->nb;
  if (EnergyEquation) {
    ComputeSoundSpeed (Rho, Energy, sys);
    /* it is necessary to update calculation of soundspeed if one uses
       alphaviscosity in FViscosity. It is not necessary in locally
       isothermal runs since sound-speed is constant. It is computed
       here for the needs of ConditionCFL. */
  }
  /* No dust particles: timestep constrained by CFL */
  if ( (NBPART == 0) || (NoTimestepConstraintByParticles) ) {
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label, DRho, DVrad, DVtheta);
      if (SloppyCFL == YES) // not the case by default
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DVrad, DVtheta, Rho, Energy, DT-dtemp);
    }
    MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    /* dthydro is the hydrodynamic timestep */
    dthydro = DT / (real)GasTimeStepsCFL;
    if (NbPlanets > 1) {
      /* dtnbody is the n-body timestep */
      dtnbody = 1e5;
      for (k = 0; k < NbPlanets; k++) {
	      for (j = k+1; j < NbPlanets; j++) {
          xk = sys->x[k];
          xj = sys->x[j];
          yk = sys->y[k];
          yj = sys->y[j];
          mk = FinalPlanetMass[k];
          mj = FinalPlanetMass[j];
          dist = sqrt( (xk-xj)*(xk-xj) + (yk-yj)*(yk-yj) );
          buf = 2.0*M_PI*sqrt( dist*dist*dist / (mk+mj+1e-8) )/BINARYSECURITY;
          dtnbody = min2(dtnbody,buf);
	      }
      }
    } else {
      /* if there is only one planet, we set dtnbody equal to dthydro */
      dtnbody = dthydro;
    }
    /* dt is the minimum between dthydro and dtnbody */
    if (IsDisk == YES) {
      dt = min2(dthydro,dtnbody);
    } else {
      dt = dtnbody;
    }
  } 
  else {
    // Case with dust particles: timestep constrained by global
    //   variable DUSTTIMESTEP entred in .par file
    if (IsDisk == YES)
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label, DRho, DVrad, DVtheta);
    dt = DustTimeStep;
    MINDT = DustTimeStep;
    MAXDT = DustTimeStep;
  }
  
  while (dtemp < 0.999999999*DT) {

    // CB: March 2022 only apply MassTaper in calculation of disc
    // potential in PframeForce.c, like in Fargo3D
    if (MASSTAPER > 1e-2) {
      MassTaper = (PhysicalTime-PhysicalTimeInitial)/(MASSTAPER*2.0*M_PI);
      MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
      if (MassTaper < 1.0) {
        for (k = 0; k < NbPlanets; k++)
          sys->mass[k] = InitialPlanetMass[k] + (FinalPlanetMass[k]-InitialPlanetMass[k])*MassTaper;
        }
      Particles_Mass = Particles_Mass_Initial*MassTaper;
    } 
    else {
      for (k = 0; k < NbPlanets; k++) {
        // if no accretion, otherwise let the planet grow!
        if (sys->acc[k] < 1e-10)
          sys->mass[k] = FinalPlanetMass[k];
      }
      Particles_Mass = Particles_Mass_Initial;
    }

    /* The current planet masses at t=PhysicalTime */
    if (NbPlanets > 1) {
      dtnbody = 1e5;
      for (k = 0; k < NbPlanets; k++) {
	      for (j = k+1; j < NbPlanets; j++) {
          xk = sys->x[k];
          xj = sys->x[j];
          yk = sys->y[k];
          yj = sys->y[j];
          mk = sys->mass[k];
          mj = sys->mass[j];
          dist = sqrt( (xk-xj)*(xk-xj) + (yk-yj)*(yk-yj) );
          buf = 2.0*M_PI*sqrt( dist*dist*dist / (mk+mj+1e-8) )/BINARYSECURITY;
          /* dtnbody is the n-body timestep */
          dtnbody = min2(dtnbody,buf);
	      }
      }
    } else 
      /* if there is only one planet, we set dtnbody to something
	 arbitrarilly large */
      dtnbody = 100.0;

    if ( (NBPART == 0) || (NoTimestepConstraintByParticles) ) {
      if (IsDisk == YES) {
	      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label, DRho, DVrad, DVtheta);
        if (SloppyCFL == NO) {  // case by default
          gastimestepcfl = 1;
          gastimestepcfl = ConditionCFL (Vrad, Vtheta, DVrad, DVtheta, Rho, Energy, DT-dtemp);
          MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
          /* dthydro is the hydrodynamic timestep */
          dthydro = (DT-dtemp)/(real)GasTimeStepsCFL;
          /* dt is the minimum between dthydro and dtnbody */
          dt = min2(dthydro,dtnbody);
        }
	      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
      } else
	      dt = dtnbody;
    } 
    else {
      if (IsDisk == YES)
        CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label, DRho, DVrad, DVtheta);
        dt = DustTimeStep;
        MINDT = DustTimeStep;
        MAXDT = DustTimeStep;
        AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    dtemp += dt;
    if (dt < MINDT)
      MINDT = dt;
    if (dt > MAXDT)
      MAXDT = dt;

    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES) GetPsysInfo (sys, MARK);

    if (IsDisk == YES) {
      /* Indirect term of star potential */
      //DiskOnPrimaryAcceleration   = ComputeAccel (force, Rho, 0.0, 0.0, 0.0, 0.0, sys);
      // CB: new nov 2025 (Sergei's proposition)
      DiskOnPrimaryAcceleration   = ComputeAccel (force, Rho, 0.0, 0.0, INDIRECTTERMSMOOTHING, 0.0, sys);
      /* Gravitational potential from star and planet(s) */
      FillForcesArrays (sys);
      /* Planets' velocities are updated with gravitationnal
	    interaction with disk */
      AdvanceSystemFromDisk (force, Rho, Energy, sys, dt);
    }

    /* Planets' positions and velocities are updated with
    gravitational interaction with star and other planets */
    AdvanceSystemRK5 (sys, dt);

    /* Dust particles' positions and velocities are updated with
    gravitational interaction with star, planets, and gas drag.
    Tests have shown that doing a simple first-order integrator
    gives in practice the same results as a second-order Leapfrog
    integrator, so in the following we first update positions with
    dt (semi-update with 2xdt timestep) and then update velocities
    with dt as well */
    if (NBPART != 0) {
      SemiUpdateDustPositions(dsys, sys, Rho, 2.*dt);  // update particles positions over dt
      UpdateDustVelocities(dsys, sys, Vrad, Vtheta, Rho, dustpcdens, dt);  // update particles velocities over dt
      /* Leapfrog integrator 
      SemiUpdateDustPositions(dsys, sys, Rho, dt);  // update particles positions over dt/2
      UpdateDustVelocities(dsys, sys, Vrad, Vtheta, Rho, dustpcdens, dt);  // update particles velocities over dt
      SemiUpdateDustPositions(dsys, sys, Rho, dt);  // update particles positions over dt/2
      */
    }

    /* Below we correct vtheta, the planet's position and velocities
       if we work in a frame non-centered on the primary */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES) {
        CorrectVtheta (Vtheta, domega);
        CorrectVtheta (DVtheta, domega);
      }
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    if (NBPART != 0) {
      RotateDsys (dsys, OmegaFrame*dt);
      /* CB 02/2024 needed to later damp pcdensX.dat files */
      if (Write_DustDensity == YES)
	      interpolation(dsys, Vrad, Vtheta, Rho, dustpcdens, dt);
    }
      
    /* Now we update gas fields */
    if (IsDisk == YES) {
      /* NEW (May 2016): evanescent boundary condition can now damp toward
      1D viscously evolving density and radial velocity profiles. Both
      profiles are solved through an independent 1D viscous evolution
      problem without planets. Nov 2025: can well be calculated even if 
      evanescent boundaries aren't used! */
      if (DampToViscous)
	      SolveOneDViscousEvolution (dt);
      
      //ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, dt, sys);

      Crashed = DetectCrash (Rho);    /* test for negative density values */
      if (DustFluid) 
        Crashed = DetectCrash (DRho);    /* test for negative density values */
      Crashed = DetectCrash (Energy);  /* test for negative energy values */
      if (Crashed == YES) {
        if (AlreadyCrashed == 0) {
          timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
          fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
          WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
          WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
          WriteDiskPolar (Vtheta, 999); /* of what happened */
          if (DustFluid) {
            WriteDiskPolar (DRho, 999);    /* We write the HD arrays */
            WriteDiskPolar (DVrad, 999);   /* in order to keep a track */
            WriteDiskPolar (DVtheta, 999); /* of what happened */
          }
          WriteDiskPolar (Energy, 999);
        }
        AlreadyCrashed++;
        masterprint ("c");
      } 
      else {
        masterprint (".");
      }
      fflush (stdout);

      if (ZMPlus) {
	      /* To model the non-axisymmetric component of the gas
	      self-gravity with an anisotropic pressure (see aniso.c) */
	      compute_anisotropic_pressurecoeff (sys);
      }

      /* Thermal diffusion needs to be applied first */
      if (EnergyEquation) {
      	if (EntropyDiffusion || RadiativeDiffusion) {
	        if (EntropyDiffusion)
	          ComputeEntropyDiffusion (Rho, Energy);
	        if (RadiativeDiffusion && !ImplicitRadiativeDiffusion)
	          ComputeRadiativeDiffusion (Rho, Energy);
	        SubStep0 (Rho, Energy, dt);
	      }
      }
      ComputeSoundSpeed (Rho, Energy, sys);
      ComputePressureField (Rho, Energy);

      if (DustFluid)
	      ComputeDustSoundSpeedAndPressure (DRho);

      /* Update vrad and vtheta with pressure, gravity and curvature
	    source terms */
      SubStep1 (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, sys, dt);

      // Modified sept 2025
      ActualiseGas (Vrad, VradInt);
      ActualiseGas (Vtheta, VthetaInt);
      if (DustFluid && (!ShortFrictionTimeApproximation)) {
        ActualiseGas (DVrad, DVradInt);
        ActualiseGas (DVtheta, DVthetaInt);
      }
      //
      //ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, dt, sys);
      //
      ActualiseGas (VradInt, Vrad);
      ActualiseGas (VthetaInt, Vtheta);
      if (DustFluid && (!ShortFrictionTimeApproximation)) {
        ActualiseGas (DVradInt, DVrad);
        ActualiseGas (DVthetaInt, DVtheta);
      }

      /* Add artifical viscosity */
      SubStep2 (Rho, Energy, DRho, dt);

      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      if (EnergyEquation)
        ActualiseGas (Energy, EnergyInt);

      if (DustFluid && (!ShortFrictionTimeApproximation)) {
        ActualiseGas (DVrad, DVradNew);    /* actualize dust vrad and vtheta */
        ActualiseGas (DVtheta, DVthetaNew);
      }

      if (DustFluid && ShortFrictionTimeApproximation) {
	      ActualiseGas (DVrad, DVradInt);    /* actualize dust vrad and vtheta */
        ActualiseGas (DVtheta, DVthetaInt);
      }

      /* Compute velocity divergence first after velocity has been updated */
      // SHOULD COME inside IF (ENERGY EQUATION) below??
      ComputeDivergenceVelocity (Vrad, Vtheta, DVrad, DVtheta);

      if (EnergyEquation) {
	      // CB: March 2022: added (Frederic's advise)
	      //ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, dt, sys);
        ActualiseGas (EnergyInt, Energy);
	
        /* Update thermal energy with heating, cooling source terms */
        if (ViscousHeating) {
          ComputeViscousTerms (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho);
          ComputeViscousHeating (Rho);
        }
        
        /* call to substep 3*/
        SubStep3 (Rho, Vtheta, dt);
        ActualiseGas (Energy, EnergyNew);
      }
      //ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, dt, sys);

      /* Update velocities, surface density and thermal energy with
      advective terms 
      CB (04/2020): with dust feedback via particles, we need first 
      to communicate boundaries!
      */
      if (DustFeedback && (IsDisk == YES))
	      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label, DRho, DVrad, DVtheta);

      Transport (Rho, Vrad, Vtheta, Energy, Label, dt);

      /* Add density and temperature floors */
      if (AddFloors) {
        AddFloorDensity (Rho);
        AddFloorEnergy (Energy);
      }

      /* Transport and diffusion when dust is modelled as a low-pressure fluid */
      if (DustFluid) {
        //if (ShortFrictionTimeApproximation) 
        //  SFTAvelocity (Rho, Vrad, Vtheta, DRho, DVrad, DVtheta);
        Transportd (DRho, Rho, DVrad, DVtheta, Label, dt);
        if (DustDiffusion == YES) 
          Diffd (DRho, Rho, DVrad, DVtheta, dt);
        if (AddFloors)
          AddFloorDensity (DRho);
      }

      /* Call to routine that reestablishes initial surface density on
	    fixed timescale */
      if (AddMass)
        DampDensity(Vrad, Vtheta, Rho, Energy, dt, sys);

      /* Simple treatment of disc evaporation by slowly decreasing the 
	    axisymmetric surface density profile of the gas */
      if (DiscEvaporation)
        Evaporation(Rho, dt);

      /* Photoevaporation via X-Ray star luminosity */
      if (PhotoEvaporation)
        ApplyPhotoEvaporation (Vrad, Rho, dt);

      /* CB (new Nov 2025): call to boundary conditions only once at the end of timestep */
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, dt, sys);

       /* Further apply so-called Stockholm damping in wave-killing zones 
       near the grid's inner and outer edges to avoid wave reflection */
      if (Evanescent)
        EvanescentBoundary (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, dt);

      /* Update of gas temperature for output */
      ComputeTemperatureField (Rho, Energy);

      /* Calculate mass (and mass excess) inside the planet's
	    circumplanetary disk */
      if (ComputeCPDMass) {
        mdcp = CircumPlanetaryMass (Rho, sys);
        exces_mdcp = mdcp - mdcp0;
      }

    } // end if disk is Yes

    PhysicalTime += dt;

  }
  masterprint ("\n");

}


void SubStep0 (Rho, Energy, dt)
     PolarGrid *Rho, *Energy;
     real dt;
{
  int i, j, l, nr, ns;
  real *energy, *entropydiff, *raddiff, *dens, *temp;
  real omegaSOR_best = omegaSOR;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  temp = Temperature->Field;
  energy = Energy->Field;
  entropydiff = EntropyDiff->Field;
  raddiff  = RadiativeDiff->Field;

  /* Update (explicite) with thermal diffusion */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (EntropyDiffusion)
	      energy[l] += entropydiff[l]*dt;
      if (RadiativeDiffusion && !ImplicitRadiativeDiffusion)
	      energy[l] += raddiff[l]*dt;
    }
  }

  /* Implicit update with radiative diffusion */
  if (ImplicitRadiativeDiffusion) {
    //niterSOR = GLOBAL_SORsolver_RadiativeDiffusion (Rho, Energy, dt, omegaSOR);
    niterSOR = SORsolver_RadiativeDiffusion (Rho, Energy, dt, omegaSOR);
    if (niterSOR < niterbuf)
      omegaSOR_best = omegaSOR;
    else {
      omegaSOR *= 0.9999;
      if ( (omegaSOR >= 2.0) || (omegaSOR <= 1.0) )
	      omegaSOR = omegaSOR_best;	
    }
    //masterprint ("niterSOR = %d, niterbuf = %d, omegaSOR = %lg, omegaBEST = %lg\n",niterSOR,niterbuf,omegaSOR,omegaSOR_best);
    niterbuf = niterSOR;
  }
}


void SubStep1 (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, sys, dt)
     PolarGrid *Vrad, *Vtheta, *Rho;
     PolarGrid *DVrad, *DVtheta, *DRho;
     PlanetarySystem *sys;
     real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  boolean selfgravityupdate;
  real *vrad, *vtheta, *rho;
  real *dvrad, *dvtheta, *drho;
  real *Pot, *Press, *cs;
  real *vradint, *vthetaint;
  real *dvradint, *dvthetaint, *dPress, *St;
  real *radfbacc, *azifbacc;
  real *indpot, *radindacc, *aziindacc;
  real *radgradp, *azigradp;
  real invomega, ts, gradp, gradphi, vt2, dvt2, dxtheta, dgradp, dforcer, dforcet;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  real Cdrag, Re, k_D, f_D, Kn, lambda, densg_cgs, Mach;
  real dust_internal_density, dust_radius, rho_code_tocgs;
  real vr_over_cs;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
  Press = Pressure->Field;
  indpot    = IndPotential->Field;
  radindacc = RadIndAcc->Field;
  aziindacc = AziIndAcc->Field;
  radfbacc = RadFBAcc->Field;
  azifbacc = AziFBAcc->Field;
  radgradp = RadGradP->Field;
  azigradp = AziGradP->Field;
  cs = SoundSpeed->Field;

  if (DustFluid) {
    drho = DRho->Field;
    dvrad = DVrad->Field; 
    dvtheta = DVtheta->Field;
    dvradint = DVradInt->Field;
    dvthetaint = DVthetaInt->Field;
    dPress = DPressure->Field;  // dust pressure calculated just before call to substep1
    St = Stokes->Field;

    /* Convert particle's mean density from g/cm^3 to code units */
    dust_internal_density = RHOPART*1000.0*pow(unit_length, 3.)/unit_mass;
    
    /* Convert particle's radius from meters to code units */
    dust_radius = Sizepart/unit_length;
    rho_code_tocgs = 0.1*unit_mass*pow(unit_length,-2.0);

    for (i = 0; i < nr; i++){
      for (j = 0; j < ns; j++){
        l = j+i*ns;

        /* Molecular mean-free path ~ m_H2 H / Sigma_gas d^2 with m_H2
          mass of H2, H pressure scale height and d ~ typical
          diameter of a particle */
        densg_cgs = rho[l]*rho_code_tocgs;  // in g cm^-2
        lambda = 3.34e-8 * pow(densg_cgs,-1.) * AspectRatio(Rmed[i]) * Rmed[i]; // in code units
        
        /* Knudsen number, note that dust_radius is already in code units */
        Kn = 0.5*lambda/dust_radius;

        /* Gas relative Mach number at particle's position: |dV| / cs */
        Mach = sqrt( (dvrad[l]-vrad[l])*(dvrad[l]-vrad[l]) + (dvtheta[l]-vtheta[l])*(dvtheta[l]-vtheta[l]) ) / cs[l];

        /* f_D coefficient to link the subsonic and the supersonic regimes */
        f_D = sqrt( 1.0 + 9.0*M_PI*Mach*Mach/128.0 );

        /* Reynolds number */
        Re = 3.0*sqrt(PI/8.0)*Mach/Kn;

        /* k_D coefficient that describes the Stokes regime */
        if (Re <= 500.0)
          k_D = 1.0+0.15*pow(Re,0.687);
        else {
          if ( (Re > 500.0) && (Re <= 1500.0) )
            k_D = 3.96e-6*pow(Re,2.4);
          if (Re > 1500.0)
            k_D = 0.11*Re;
        }
        
        /* Final expression for drag coefficient Cdrag as in Paardekooper 07 */
        //Cdrag = pow(3.0*Kn+1.0,2.0) * pow(9.0*Kn*Kn*f_D + 3.0*Kn*k_D,-1.0);
        Cdrag = (3.0*Kn+1.0)*(3.0*Kn+1.0) / (9.0*Kn*Kn*f_D + 3.0*Kn*k_D);
        /* Dimensionless stopping time at the particle = stokes number */
        St[l] =  0.5*M_PI*Cdrag*dust_internal_density*dust_radius/rho[l];
      }
    }
  } 
  /* In this substep we take into account the source terms of Euler
     equations. We update velocities with pressure gradients,
     gravitational forces and curvature terms */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      invomega = pow(Rmed[i],1.5);
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        gradp = (Press[l]-Press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
        radgradp[l] = (Press[l]-Press[lim])*InvDiffRmed[i];
        gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
        radindacc[l] = -(indpot[l]-indpot[lim])*InvDiffRmed[i];
        vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
        vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
        vt2 = vt2*vt2;
        vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
        if (DustFeedback && !DustFluid) {
          vradint[l] += dt*radfbacc[l];
        }
        if (DustFluid) {
          ts = St[l] * invomega;
          if (!ShortFrictionTimeApproximation) {
            dgradp = (dPress[l]-dPress[lim])*2.0/(drho[l]+drho[lim])*InvDiffRmed[i];
            dvt2 = dvtheta[l]+dvtheta[ljp]+dvtheta[lim]+dvtheta[ljp-ns];
            dvt2 = dvt2/4.0+Rinf[i]*OmegaFrame;
            dvt2 = dvt2*dvt2;
            dvradint[l] = dvrad[l]+dt*(-gradphi-dgradp+dvt2*InvRinf[i]);
            if (DustFeelDisk) {
              //dforcer   = -(dvrad[l]-vrad[l])/(St[l]*pow(Rmed[i],1.5) + dt/2.);
              //dforcer = -(dvrad[l]-vrad[l])/(St[l]*pow(Rmed[i],1.5));
              //dvradint[l] += dt*dforcer;
              // implicite scheme test:
              dvradint[l] = (dvradint[l]*ts + vrad[l]*dt)/(ts+dt);
            }
          }
          /* Dust velocity in the short-friction time approximation: original
            expression from Johansen & Klahr (2005), generalised to include
            dust drag on gas */
          else {
            if (DustFeelDisk) {
              if (DustFeedback)
          ts /= (1.0+drho[l]/rho[l]);
              dvradint[l] = vrad[l] + ts*gradp;
            }
          }
          if (DustFeedback) {
            ts = St[l] * invomega;
            //radfbacc[l] = (drho[l]/rho[l])*(dvrad[l]-vrad[l])/(St[l]*pow(Rmed[i],1.5));
            radfbacc[l]   = (drho[l]/rho[l])*(dvrad[l]-vrad[l])/(ts + dt/2.);
            //vradint[l] += dt*radfbacc[l];
            // implicite scheme test:
            vradint[l] = (vradint[l]*ts + (drho[l]/rho[l])*dvrad[l]*dt) / (ts + dt*(drho[l]/rho[l]));
          }
        }
      }
    }

#pragma omp for
    for (i = 0; i < nr; i++) {

      dxtheta = (PMAX-PMIN)/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      
      //supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      // CB (May 2019): changed expression of supp_torque so that it applies to an 
      // arbitrary initial density profile (not necesarilly a power-law function or R)
      invomega = pow(Rmed[i],1.5);
      supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5)*SIGMA0/SigmaMed[i];
      if (CavityTorque) {
        /* June 2022 (Guillaume Robert's internship) */
        vr_over_cs = -( 0.5*(FRACINT+FRACEXT) + 0.5*(FRACINT-FRACEXT)*tanh((CAVITYRADIUS-Rmed[i])/CAVITYWIDTH) );  // negative definite
      }
      
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        gradp = (Press[l]-Press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
        azigradp[l] = (Press[l]-Press[ljm])*invdxtheta;
        if ( ZMPlus ) {
          /* To model the non-axisymmetric component of the gas
            self-gravity with an anisotropic pressure (see aniso.c) */
          gradp *= SG_aniso_coeff;
        }
        gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
        aziindacc[l] = -(indpot[l]-indpot[ljm])*invdxtheta;
        vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
        /* June 2022 (Guillaume Robert's internship) */
        if (CavityTorque)
          supp_torque = 0.5*pow(Rmed[i],-1.5)*vr_over_cs*cs[l];
        /* --------- */
        vthetaint[l] += dt*supp_torque;
        if (DustFeedback && !DustFluid) {
          vthetaint[l] += dt*azifbacc[l];
        }
        if (DustFluid) { 
          ts = St[l] * invomega;
          if (!ShortFrictionTimeApproximation) {
            dgradp = (dPress[l]-dPress[ljm])*2.0/(drho[l]+drho[ljm])*invdxtheta;
            dvthetaint[l] = dvtheta[l]-dt*(gradphi+dgradp);
            if(DustFeelDisk){
              //dforcet =   -(dvtheta[l]-vtheta[l])/(St[l]*pow(Rmed[i],1.5) + dt/2.);
              //dforcet = -(dvtheta[l]-vtheta[l])/(St[l]*pow(Rmed[i],1.5));
              //dvthetaint[l] += dt*dforcet; 
              // implicite scheme test:
              dvthetaint[l] = (dvthetaint[l]*ts + vtheta[l]*dt)/(ts+dt);
            }
          }
          /* Dust velocity in the short-friction time approximation: original
            expression from Johansen & Klahr (2005), generalised to include
            dust drag on gas */
          else {
            if (DustFeelDisk) {
              if (DustFeedback)
          ts /= (1.0+drho[l]/rho[l]);
              dvthetaint[l] = vtheta[l] + ts*gradp;
            }
          }
          if (DustFeedback) {
            ts = St[l] * invomega;
            //azifbacc[l] = (drho[l]/rho[l])*(dvtheta[l]-vtheta[l])/(St[l]*pow(Rmed[i],1.5));
            azifbacc[l] =   (drho[l]/rho[l])*(dvtheta[l]-vtheta[l])/(ts + dt/2.);
            //vthetaint[l] += dt*azifbacc[l];
            // implicite scheme test:
            vthetaint[l] = (vthetaint[l]*ts + (drho[l]/rho[l])*dvtheta[l]*dt)/(ts+dt*(drho[l]/rho[l]));
          }
        }
      }
    }
  }
  
  /* Update velocities with gas self-gravity */
  if ( SelfGravity ) {
    selfgravityupdate = YES;
    compute_selfgravity(Rho, VradInt, VthetaInt, DVradInt, DVthetaInt, dt, selfgravityupdate);
  }
  
  /* Update velocities with gas viscosity */
  if ( (ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0) || ( DustFluid && ( (DALPHAVISCOSITY != 0.0) || (DVISCOSITY != 0.0)) ) ) {
    ComputeViscousTerms (VradInt, VthetaInt, Rho, DVradInt, DVthetaInt, DRho);
    UpdateVelocitiesWithViscosity (VradInt, VthetaInt, Rho, DVradInt, DVthetaInt, DRho, dt);
  }
  
}


void SubStep2 (Rho, Energy, DRho, dt)
     PolarGrid *Rho, *Energy, *DRho;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy, *energyint;
  real *dvrad, *dvtheta, *drho;
  real *vradnew, *vthetanew, *dvradnew, *dvthetanew, *qt, *qr, *dqt, *dqr, *test;
  real dxtheta, invdxtheta;
  real dv, ddv;
  
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;
  test = Test->Field;
  //
  if (DustFluid) {
    drho = DRho->Field;
    dvrad = DVradInt->Field;
    dvtheta = DVthetaInt->Field;
    dqr = DRhoInt->Field;
    dvradnew = DVradNew->Field;
    dvthetanew = DVthetaNew->Field;
    dqt = DTemperInt->Field;
  }
  //
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0;
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	      qt[l] = 0.0;
      /* Dust Fluid */
      if (DustFluid) {
	      ddv = dvrad[lip]-dvrad[l];
        if (ddv < 0.0)
          dqr[l] = CVNR*CVNR*drho[l]*ddv*ddv;
        else
          dqr[l] = 0.0;
        ddv = dvtheta[ljp]-dvtheta[l];
        if (ddv < 0.0) 
          dqt[l] = CVNR*CVNR*drho[l]*ddv*ddv;
        else
          dqt[l] = 0.0;
      }
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
        if (DustFluid)
          dvradnew[l] = dvrad[l]-dt*2.0/(drho[l]+drho[lim])*(dqr[l]-dqr[lim])*InvDiffRmed[i]; /* dust */
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = (PMAX-PMIN)/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lip = l+ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
        if (DustFluid)
          dvthetanew[l] = dvtheta[l]-dt*2.0/(drho[l]+drho[ljm])*(dqt[l]-dqt[ljm])*invdxtheta; /* dust */
        if (EnergyEquation) {
          energyint[l] = energy[l] -				\
            dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -		\
            dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
          /*
          test[l] = -qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -	\
            qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
          */
        }
      }
    }
  }
}

	       
void SubStep3 (Rho, Vtheta, dt)
     PolarGrid *Rho, *Vtheta;
     real dt;
{
  extern boolean TempPresc, BetaCooling, StellarIrradiation, ImposedStellarLuminosity;
  int i, j, l, nr, ns;
  real *energy, *energynew, *dens, *divergence, *vischeat, *thercool, *vtheta, *fbedot, *opacity;
  real num, den, omega, beta, coolingtime, tau, tau_eff, temp, temp_irr;
  real dTdRp, dTdRm, dTdR, costheta, buf;
  real albedo=0.1, Lstar;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  vtheta = Vtheta->Field;
  energy = EnergyInt->Field;
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  vischeat = ViscHeat->Field;
  thercool = ThermCool->Field;
  fbedot   = FBedot->Field;
  opacity = Opacity->Field;

  /* In this substep, we update the gas thermal energy with source
     terms (compression/dilatation, viscous heating, radiative
     (thermal) cooling, simple beta cooling, or via a temperature
     prescription) */
#pragma omp parallel private(j,l,num,den)
  {
#pragma omp for

    if (!ThermalCooling) {
      /* First Update (implicite) with pdV only */
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          num = energy[l] + (ViscousHeating == YES ? 1 : 0)*vischeat[l]*dt;
          den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l];
          energynew[l] = num/den;
	      }
      }
      /* Here cooling is modeled through a temperature prescription,
	    with characteristic time is PrescTimeMed(r). Energy update is
	    implicite. */
      if (TempPresc) {
        for (i = 0; i < nr; i++) {
          for (j = 0; j < ns; j++) {
            l = j+i*ns;
            num = energynew[l] + EnergyMed[i]*(dens[l]/SigmaMed[i])*(dt/PrescTimeMed[i]);
            den = 1.0+dt/PrescTimeMed[i]; 
            energynew[l] = num/den;
          }
        }
      }
      /* Case of a simple beta-cooling, where the cooling source term
	     reads -e/tau, with tau = beta/Omega */
      if (BetaCooling) {
        for (i = 0; i < nr; i++) {
          omega = 0.0;
          for (j = 0; j < ns; j++) {
            l = j+i*ns;
            omega += (vtheta[l] + Rmed[i]*OmegaFrame);
          }
          omega /= (real)ns;
          omega /= Rmed[i];
          beta  = ComputeBetaCooling(Rmed[i]);
          for (j = 0; j < ns; j++) {
            l = j+i*ns;
            num = energynew[l];
            coolingtime = beta / omega;
            den = 1.0+dt/coolingtime; 
            energynew[l] = num/den;
          }
        }
      }
    }

    else { // if thermal cooling set to yes
      /* Implicite update with thermal cooling, viscous heating and
	    stellar irradiation in option */
      if (StellarIrradiation && ImposedStellarLuminosity) {
        // convert stellar luminosity from Watt to code units
        // remember that STELLARLUMINOSITY is set in Solar luminosities in .par file
        Lstar = STELLARLUMINOSITY*3.8e26 * pow(unit_mass,-1.0) * pow(unit_length,-2.0) * pow(unit_time, 3.0);
      }
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          tau = 0.5*opacity[l]*dens[l];      // calculation valid even if a constant opacity is set in .par file
          tau_eff = 0.375*tau + 0.25*sqrt(3.0) + 0.25/tau; // effective optical depth
          temp = (ADIABATICINDEX-1.0)*energy[l]*pow(dens[l],-1.0);  // temperature at time t

          if (StellarIrradiation) {

            if (!ImposedStellarLuminosity) {
              /* case where stellar irradiation sets a prescribed
              background temperature profile, whose value at code
              unit of length and power-law exponent are set in the
              .par parameter file */
              // BACKGROUNDTEMPERATURE is meant to be the gas temperature 
              // set by stellar irradiation at code's unit of length
              temp_irr = (BACKGROUNDTEMPERATURE/unit_temperature) * pow(Rmed[i],SLOPEBACKGROUNDTEMPERATURE);
            } else {
              /* case where irradiation heating term calculated with
              imposed stellar luminosity. In that case the heating
              source term reads (1-beta) Lstar/4piR^2 cos(theta)
              with beta the albedo (coded in raw in current
              routine), Lstar the stellar luminosity set in Solar
              luminosities in the .par parameter file, and theta
              the grazing angle. We work out the background
              temperature profile temp_irr that is equivalent to
              above source term */
              // (i) cosinus of grazing angle theta: cos(theta) = dH/dR - H/R
              // with H = T^1/2 / Omega. This gives: 
              // cos(theta) = 0.5sqrt(RT(R)) x (1 + dlog(T)/dlog(R))
              if ( (i > 0) && (i < nr-1) ) {
                dTdRp = InvDiffRmed[i+1]*((ADIABATICINDEX-1.0)*energy[l+ns]*pow(dens[l+ns],-1.0)-temp);
                dTdRm = InvDiffRmed[i]*(temp-(ADIABATICINDEX-1.0)*energy[l-ns]*pow(dens[l-ns],-1.0));
                dTdR  = 0.5*(dTdRp+dTdRm);
              } else {
              if (i==0)
                dTdR = InvDiffRmed[i+1]*((ADIABATICINDEX-1.0)*energy[l+ns]*pow(dens[l+ns],-1.0)-temp);
              if (i==nr-1)
                dTdR = InvDiffRmed[i]*(temp-(ADIABATICINDEX-1.0)*energy[l-ns]*pow(dens[l-ns],-1.0));
              }
              costheta = 0.5*sqrt(Rmed[i]*temp)*(1.0 + Rmed[i]*dTdR/temp);
              if (costheta > 1) {
                //printf ("Issue with costheta calculation in substep3 as costheta = %lg\n",costheta);
                //printf ("i = %d, tau_eff = %lg, shadow = %lg, costheta = %lg, temp = %lg, dTdT = %lg\n",i, tau_eff, shadow_function[j], costheta, temp, dTdR);
                //prs_exit();
                costheta = 1.0;
              }
              if (costheta < 0) {
                //printf ("Issue with costheta calculation in substep3 as costheta = %lg\n",costheta);
                //printf ("i = %d, tau_eff = %lg, shadow = %lg, costheta = %lg, temp = %lg, dTdT = %lg\n",i, tau_eff, shadow_function[j], costheta, temp, dTdR);
                costheta = 0.0;
              }
              //costheta = 0.04;  // cuidadin!
              costheta = 0.5*sqrt(Rmed[i]*temp)*(2.0*FLARINGINDEX);  // cuidadin!
              // (iii) obtain equivalent irradiation temperature profile:
              buf = 0.5*(tau_eff/sigma_SB)*(1.0-albedo)*(Lstar/4.0/M_PI/Rmed[i]/Rmed[i])*costheta;
              temp_irr = pow(buf,0.25);
              //
              //printf("i = %d: tau_eff = %lg, costheta = %lg, buf = %lg, temp_irr = %lg, temp = %lg\n", i, tau_eff, costheta, buf, temp_irr, temp);
              //if ((PhysicalTime != 0.0) && (i == nr-1))
              //	prs_exit();
              //
              if (isnan(temp_irr) == 1) {
                masterprint ("NaN issue with temp_irr calculation in substep 3\n");
                masterprint ("buf = %lg, i = %d, tau_eff = %lg, costheta = %lg, temp = %lg\n",buf, i, tau_eff, costheta, temp*unit_temperature);
                prs_exit();
              }
            }
            
            num = energy[l] + dt*(vischeat[l]*(ViscousHeating == YES ? 1 : 0) +	\
                (2.0*sigma_SB/tau_eff)*(3.0*pow(temp,4.0) + pow(temp_irr,4.0)));
          } 
          else {  // If no stellar irradiation
            num = energy[l] + dt*(vischeat[l]*(ViscousHeating == YES ? 1 : 0) +	\
                (6.0*sigma_SB/tau_eff)*pow(temp,4.0));
          }
          den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l]+		\
            8.0*dt*(sigma_SB/tau_eff)*(ADIABATICINDEX-1.0)*pow(temp,3.0)*pow(dens[l],-1.0);
          energynew[l] = num/den;
        }
      }
    }
    /* Update (explicite) with heating due to dust feedback on gas
       (beta) */
    if (DustFeedback) {
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          energynew[l] += fbedot[l]*dt;
	      }
      }
    }
  }
}
 
int ConditionCFL (Vrad, Vtheta, DVrad, DVtheta, Rho, Energy, deltaT)
     PolarGrid *Vrad, *Vtheta;
     PolarGrid *DVrad, *DVtheta;
     PolarGrid *Rho, *Energy;
     real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  static real DVresidual[MAX1D], DVmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  int ideb=0, jdeb=0;
  real invdt1, invdt2, invdt3, invdt4, invdt5, invdt6, cs, newdt, dt, dtg, dtd, dt_fb, dtshear, dt1D, buf;
  real invdt7, invdt8, invdt9, invdt10, invdt11, invdt12;
  real itdbg1, itdbg2, itdbg3, itdbg4, itdbg5, itdbg6, mdtdbg; /* debugging variables */
  real itdbg7, itdbg8, itdbg9, itdbg10, itdbg11, itdbg12; /* debugging variables */
  real dphi, invdphi, gradT, Rl, lambda, chicoeff;
  real viscosity=0.0, dviscosity=0.0;
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr, visct;
  real *dustvr, *dustvt, *St, csd, ddvr, ddvt;
  real *soundspeed, *dsoundspeed, *temperature, *dens, *opacity;
  extern boolean Evanescent, DampToViscous;

  if (EnergyEquation && RadiativeDiffusion && !ImplicitRadiativeDiffusion)
     ComputeOpacities (Rho, Energy);

  opacity = Opacity->Field;
  soundspeed = SoundSpeed->Field;
  temperature = Temperature->Field;
  dens = Rho->Field;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;

  if (DustFluid) {
    dustvt = DVtheta->Field;
    dustvr = DVrad->Field;
    St = Stokes->Field;
    dsoundspeed = DSoundSpeed->Field;
  }
  newdt = 1e30;
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*(PMAX-PMIN)/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
    //
    if (DustFluid) {
      DVmoy[i] = 0.0;
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        DVmoy[i] += dustvt[l];
      }
      DVmoy[i] /= (real)ns;
    }
    //
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    if ((ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0))
      viscosity = FViscosity(Rmed[i]);
    if ( DustFluid && ( (DALPHAVISCOSITY != 0.0) || (DVISCOSITY != 0.0)) )
      dviscosity = DFViscosity(Rmed[i]);
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*(PMAX-PMIN)/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
	      Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
	      Vresidual[j] = vt[l];	       /* Standard algorithm */
       //
      if (DustFluid) {
        if (FastTransport == YES)
          DVresidual[j] = dustvt[l]-DVmoy[i];  /* FARGO algorithm */
        else
          DVresidual[j] = dustvt[l];	       /* Standard algorithm */
      }
    }
    Vresidual[ns]=Vresidual[0];
    if (DustFluid)
      DVresidual[ns]=DVresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = soundspeed[l];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      if (viscosity != 0.0) 
	      invdt5 = viscosity*4.0/pow(min2(dxrad,dxtheta),2.0);
      else 
	      invdt5 = 1e-10;
      if (EnergyEquation && (EntropyDiffusion || (RadiativeDiffusion && !ImplicitRadiativeDiffusion))) {
        if (EntropyDiffusion)
          invdt6 = DIFFUSIVITY*4.0/pow(min2(dxrad,dxtheta),2.0);
        if (RadiativeDiffusion && !ImplicitRadiativeDiffusion) {
          gradT = sqrt(pow(InvDiffRmed[i+1]*(temperature[lip]-temperature[l]),2.0)+ pow(InvRmed[i]*invdphi*(temperature[ljp]-temperature[l]),2.0));
          Rl = (8.0*sqrt(ADIABATICINDEX*temperature[l])*pow(Rmed[i],1.5)/dens[l]/opacity[l])*gradT/temperature[l];
          if (Rl <= 2.0) {
            lambda = 2./(3.+sqrt(9.0+10.0*Rl*Rl));
          } else {
            lambda = 10./(10.*Rl+9.0+sqrt(180.*Rl+81.));
          }
          /* radiative diffusion coefficient at i,j (cell-centred) */
          chicoeff = 64.0*lambda*sigma_SB*ADIABATICINDEX*(ADIABATICINDEX-1.0)*pow(temperature[l],4.0)/opacity[l]/pow(dens[l],2.0)*pow(Rmed[i],3.0);
          invdt6 = chicoeff*4.0/pow(min2(dxrad,dxtheta),2.0);
          //invdt6 *= 2.0;  // extra precaution!
        }
      }
      else 
	      invdt6 = 1e-10;
      //
      if (DustFluid) {
        csd = dsoundspeed[l];
        if(!ShortFrictionTimeApproximation) {
          invdt7 = csd/(min2(dxrad,dxtheta));
          invdt8 = fabs(dustvr[l]+1e-15)/dxrad;
          invdt9 = fabs(DVresidual[j])/dxtheta;
          if(DustFeelDisk==NO) 
            invdt10 = 1.e-10;
          else
            invdt10 = 1./(St[l]+1.e-10)*sqrt(G*1.0/Rmed[i])/Rmed[i];
        } else {
          invdt7=1.e-10;
          invdt8=1.e-10;
          invdt9=1.e-10;
          invdt10=1.e-10;
        }
        ddvr = dustvr[lip]-dustvr[l];
        ddvt = dustvt[ljp]-dustvt[l];
        if (ddvr >= 0.0) ddvr = 1e-10;
        else ddvr = -ddvr;
        if (ddvt >= 0.0) ddvt = 1e-10;
        else ddvt = -ddvt;
        invdt11 = max2(ddvr/dxrad,ddvt/dxtheta);
        invdt11*= 4.0*CVNR*CVNR;
        if (dviscosity != 0.0) 
          invdt12 = dviscosity*4.0/pow(min2(dxrad,dxtheta),2.0);
        else 
          invdt12 = 1e-10;
      }
      //
      if (DustFluid) {
        dtg = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4+invdt5*invdt5+invdt6*invdt6);
        dtd = CFLSECURITY/sqrt(invdt7*invdt7+invdt8*invdt8+invdt9*invdt9+invdt10*invdt10+invdt11*invdt11+invdt12*invdt12);
        dt=min2(dtg,dtd);
        dt=max2(dt,1.e-10);
      } else {
	      dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4+invdt5*invdt5+invdt6*invdt6);
      }
      if (dt < newdt) {
        newdt = dt;
        ideb = i;
        jdeb = j;
        if (debug == YES) {
          ideb = i;
          jdeb = j;
          itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4; itdbg5=1.0/invdt5; itdbg6=1.0/invdt6;
          if (DustFluid)
            itdbg7 = 1.0/invdt7; itdbg8=1.0/invdt8; itdbg9=1.0/invdt9; itdbg10=1.0/invdt10; itdbg11=1.0/invdt11; itdbg12=1.0/invdt12;
          mdtdbg = newdt;
          viscr = dxrad/dvr/4.0/CVNR/CVNR;     
          visct = dxtheta/dvt/4.0/CVNR/CVNR;
        }
      }  
    }
  }
  dtshear = 1e30;
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = (PMAX-PMIN)*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < dtshear) dtshear = dt;
    if (dt < newdt) newdt = dt;
  }
  if (DustFluid && (!ShortFrictionTimeApproximation)) {
    for (i = Zero_or_active; i < MaxMO_or_active; i++) {
      dt = (PMAX-PMIN)*CFLSECURITY/(real)NSEC/fabs(DVmoy[i]*InvRmed[i]-DVmoy[i+1]*InvRmed[i+1]);
      if (dt < newdt) newdt = dt;
    }
  }
  /* NEW Jan 16: if dust feedback on gas, add dt constraint that the
     minimum stopping time of the dust particles should be resolved by
     2 hydro time steps */
  if ( (NBPART != 0) && (DustFeedback) ) {
    dt_fb = 0.5*Minimum_Stopping_Time;
    //printf ("dt_fb = %lg\n",dt_fb);
    if (dt_fb == 0.0) 
      dt_fb = 1e10;
    if (dt_fb < newdt) newdt = dt_fb;
  }
  /* NEW June 17: include dt constraint from 1D grid */
  if ( (Evanescent && DampToViscous) || DampToViscous) {
    dt1D = 1e10;
    for (i=0; i<NRAD1D; i++) {
      buf = pow(Rsup1D[i]-Rinf1D[i],2.0)/8.0/viscosity1D[i];
      if (buf < dt1D) dt1D = buf;
    }
    if (dt1D < newdt) newdt = dt1D;
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Artificial viscosity limit     : %g\n", itdbg4);
    printf ("Viscosity limit                : %g\n", itdbg5);
    printf ("Thermal or radiation Diffusivity limit      : %g\n", itdbg6);
    printf ("Radiative over viscous diffusion times      : %g\n", itdbg6/itdbg5);   
    printf ("Radiative diff over sound crossing time     : %g\n", itdbg6/itdbg1);   
    printf ("Fargo shear limit              : %g\n", dtshear);
    if (DustFluid) {
      printf ("Dust sound speed limit              : %g\n", itdbg7);
      printf ("Dust Radial motion limit            : %g\n", itdbg8);
      printf ("Dust Residual circular motion limit : %g\n", itdbg9);
      printf ("Dust Radial aerodrag limit     : %g (dvr=%g, vr=%g, St=%g)\n", itdbg10, dustvr[jdeb+ideb*ns], vr[jdeb+ideb*ns], St[jdeb+ideb*ns]);
      printf ("Dust artificial viscosity limit     : %g\n", itdbg11);
      printf ("Dust viscosity limit                : %g\n", itdbg12);
    }
    if (Evanescent && DampToViscous) {
      printf ("CFL limit due to 1D grid     : %g\n", dt1D);
    }
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (int)(ceil(deltaT/newdt));
}

void ComputeViscousHeating (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  int lip, li2p;
  real r, rip, ri2p, qpip, qpi2p, viscosity=0.0;
  real *dens, *divergence, *Trr, *Trp, *Tpp, *vischeat;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  divergence = DivergenceVelocity->Field;
  vischeat = ViscHeat->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  /* We calculate the heating source term from i=1 */
  for (i = 1; i < nr; i++) {     /* Trp defined from i=1 */
    if ((ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0))
      viscosity = FViscosity(Rmed[i]);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (viscosity != 0.0) {
        /* Note that the factor 2 in front of Trp was originally
          missing in the online version of the code (that factor 2
          was also missing in D'Angelo+ 03...) */
        vischeat[l] = 0.5/viscosity/dens[l]*( Trr[l]*Trr[l] +		\
                      2.0*Trp[l]*Trp[l] +	\
                      Tpp[l]*Tpp[l] );
        vischeat[l] += (2.0/9.0)*viscosity*dens[l]*divergence[l]*divergence[l];
      }
      else
	      vischeat[l] = 0.0;
    }
  }
  /* We calculate the heating source term Vischeat for i=0 */
  i = 0;
  r    = Rmed[i];
  rip  = Rmed[i+1];
  ri2p = Rmed[i+2];
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    lip = l+ns;
    li2p = lip+ns;
    qpip = vischeat[lip];   // vischeat(i=1,j)
    qpi2p = vischeat[li2p]; // vischeat(i=2,j)
    if (viscosity != 0.0) {
      // power-law extrapolation
      vischeat[l] = qpip*exp( log(qpip/qpi2p) * log(r/rip) / log(rip/ri2p) );
    }
    else
      vischeat[l] = 0.0;
  }
}

void ComputeOpacities (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *energ;
  real *dens, rho3D, phys_dens;
  real *opacity;
  real *test;
  real temp, phys_temp;
  real roversigma, buf;
  real temp_transition_34, temp_transition_45, temp_transition_56, temp_transition_67, temp_transition_78;
  extern boolean SetConstantOpacity;
  dens = Rho->Field;
  energ = Energy->Field;
  opacity = Opacity->Field;
  test = Test->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  if (!SetConstantOpacity) {

    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
      	l = i*ns + j;

        /* Convert code temperature into Kelvins */
        temp = (ADIABATICINDEX-1.0)*energ[l]*pow(dens[l],-1.0);  // temperature in code units
        phys_temp = temp * unit_temperature;

        /* Convert 3D volume density into g.cm^-3 */
        roversigma = Rmed[i] / dens[l];

        buf = ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]*pow(roversigma,3.);
        rho3D = 0.5*pow(buf,-0.5);  // 3D density = sigma / 2H, in code units

        phys_dens = rho3D * unit_mass * pow(unit_length, -3.);  // in kg.m^(-3)
        phys_dens *= 1e-3;  // in g.cm^(-3)
        
        /* Opacities are calculated in cm^2/g from Bell and Lin (94) tables  */
        if ( phys_temp < 167.0 )
          opacity[l] = 2e-4*pow(phys_temp,2.0);
        else {
          if ( phys_temp < 203.0 )
            opacity[l] = 2e16*pow(phys_temp,-7.0);
          else {
            temp_transition_34 = pow(2e82*phys_dens,2./49);
            if ( phys_temp < temp_transition_34 ) 
              opacity[l] = 0.1*pow(phys_temp,0.5);
            else {
              temp_transition_45 = pow(2e89*pow(phys_dens,1./3),1./27);
              if ( phys_temp < temp_transition_45 )
               opacity[l] = 2e81*pow(phys_dens,1.0)*pow(phys_temp,-24.);
              else {
                temp_transition_56 = pow(1e28*pow(phys_dens,1./3),1./7);
                if ( phys_temp < temp_transition_56 )
                  opacity[l] = 1e-8*pow(phys_dens,2./3)*pow(phys_temp,3.);
                else {
                  temp_transition_67 = pow(1.5e56*pow(phys_dens,2./3),0.08);
                  if ( phys_temp < temp_transition_67 )
                    opacity[l] = 1e-36*pow(phys_dens,1./3)*pow(phys_temp,10.);
                  else {
                    temp_transition_78 = pow(4.31e20*phys_dens,2./5);
                    if ( phys_temp < temp_transition_78 )
                      opacity[l] = 1.5e20*pow(phys_dens,1.)*pow(phys_temp,-2.5);
                    else
                      opacity[l] = 0.348;
                  }
                }
              }
            }
          }
        }
        opacity[l] *= FACTOROPACITIES;  // optional fudge factor (default = 1)
        //test[l] = opacity[l];   // in cm^2 / g
        /* We convert opacities them in m^2 / kg, before translating
          the result into code units */
        opacity[l] *= (0.1 * pow(unit_length,-2.0) * pow(unit_mass,1.0));
      }
    }
  } 

  else {
    /* Case where we impose a constant opacity throughout the grid */
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        opacity[l] = IMPOSEDCONSTANTOPACITY;  // in cm^2 / g
        /* We convert opacities them in m^2 / kg, before translating
          the result into code units */
        opacity[l] *= (0.1 * pow(unit_length,-2.0) * pow(unit_mass,1.0));
      }
    }
  }
}

/* old explicit update with thermal cooling, now changed to an
   implicit update directly in substep 3 */
/*
void ComputeThermalCooling (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *energ, *thercool, *opacity;
  real temp, tau, tau_eff;
  real temp_irr;
  extern boolean StellarIrradiation;
  ComputeOpacities (Rho, Energy);
  dens = Rho->Field;
  energ = Energy->Field;
  thercool = ThermCool->Field;
  opacity = Opacity->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      tau = 0.5*opacity[l]*dens[l]; 
      tau_eff = 0.375*tau + 0.25*sqrt(3.0) + 0.25/tau; // effective optical depth
      temp = (ADIABATICINDEX-1.0)*energ[l]*pow(dens[l],-1.0);  // temperature
      if (!StellarIrradiation)
	      thercool[l] = 2.0*sigma_SB*pow(temp,4.)*pow(tau_eff,-1.0);
      else {
        //temp_irr = pow(ASPECTRATIO,2.0)*pow(Rmed[i],-1.0+2.0*FLARINGINDEX);
        temp_irr = (BACKGROUNDTEMPERATURE/unit_temperature) * pow(Rmed[i],SLOPEBACKGROUNDTEMPERATURE);   // BACKGROUNDTEMPERATURE is meant to be the gas temperature set by stellar irradiation at code's unit of length
        thercool[l] = 2.0*sigma_SB*(pow(temp,4.)-pow(temp_irr,4.))*pow(tau_eff,-1.0);
      }
    }
  }
}
*/


void ComputeEntropyDiffusion (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  int lip, lim, ljp, ljm;
  real *energy, *tempint, *dens, *entropydiff;
  real dphi, invdphi, laplacien, buf;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* tempint is an intermediary array: energy / density */
  tempint = TempInt->Field;
  energy = Energy->Field;
  dens = Rho->Field;
  entropydiff = EntropyDiff->Field;
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  /* Thermal diffusion implemented as in Paardekooper, Baruteau & Kley
     2011 (see their Eq. 1) with actually entropy diffusion. This
     avoids the background disc structure to evolve
     significantly... */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      buf = (ADIABATICINDEX-1.0)*energy[l]/pow(dens[l],ADIABATICINDEX);
      if (buf <= 0.0) {
	      tempint[l] = 0.0;
      } 
      if (buf > 0.0) {
	      tempint[l] = log( buf );
      }
    }
  }
  for (i = 1; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      laplacien = InvDiffRsup[i]*InvRmed[i]*(
					     Rsup[i]*InvDiffRmed[i+1]*(tempint[lip]-tempint[l]) - \
					     Rinf[i]*InvDiffRmed[i]*(tempint[l]-tempint[lim]) 
					     ) +			\
	InvRmed[i]*InvRmed[i]*invdphi*invdphi*(tempint[ljp]+tempint[ljm]-2.0*tempint[l]);
      entropydiff[l] = energy[l]*DIFFUSIVITY*laplacien;
    }
  }
}

void ComputeRadiativeDiffusion (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  /* NEW (Sept. 2017): this new module solves the radiative source
term in the energy equation with an explicit scheme. In 2D, the
equation we solve is de/dt = 2H div(K grad(T)) with H the pressure
scale height and K = 16sigma_SB lambda T^3 / rho kappa, where lambda
the flux limiter as in Kley (1989), rho = Sigma/2H and kappa the
Rossland opacity (our standard opacity) */
  int i, j, l, nr, ns;
  int lip, lim, ljp, ljm;
  real *energy, *tempint, *kcoeff, *dens, *raddiff, *opacity;
  real dphi, invdphi, divF, buf, lambda, Rl, gradT;
  real KRsup, KRinf, KPhisup, KPhiinf;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* Compute opacities first, used to compute the radiative diffusion
     coefficient next */
  ComputeOpacities (Rho, Energy);
  opacity = Opacity->Field;
  tempint = TempInt->Field;
  kcoeff = RadiativeKCoeff->Field;
  energy = Energy->Field;
  dens = Rho->Field;
  raddiff = RadiativeDiff->Field;
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* temperature at i,j */
      tempint[l] = (ADIABATICINDEX-1.0)*energy[l]/dens[l];
    }
  }
  for (i = 0; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      /* flux-limiting coefficient lambda, as in Kley (1989) */
      gradT = sqrt(pow(InvDiffRmed[i+1]*(tempint[lip]-tempint[l]),2.0) + pow(InvRmed[i]*invdphi*(tempint[ljp]-tempint[l]),2.0));
      Rl = (8.0*sqrt(ADIABATICINDEX*tempint[l])*pow(Rmed[i],1.5)/dens[l]/opacity[l])*(gradT/tempint[l]);
      if (Rl <= 2.0) {
        lambda = 2./(3.+sqrt(9.0+10.0*Rl*Rl));
      } else {
        lambda = 10./(10.*Rl+9.0+sqrt(180.*Rl+81.));
      }
      /* radiative diffusion coefficient at i,j (cell-centred) */
      kcoeff[l] = 32.0*sigma_SB*lambda*sqrt(ADIABATICINDEX)*pow(tempint[l],3.5)*pow(Rmed[i],1.5)/dens[l]/opacity[l];
    }
  }
  // Special case i=nr-1 (required for calculation of div(F) below)
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    kcoeff[l] = kcoeff[l-ns];  // I could otherwise do a power-law extrapolation!
  }
  for (i = 1; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      KRsup = ((Rmed[i+1]-Rsup[i])*kcoeff[l] + (Rsup[i]-Rmed[i])*kcoeff[lip])/(Rmed[i+1]-Rmed[i]);
      KRinf = ((Rmed[i]-Rinf[i])*kcoeff[lim] + (Rinf[i]-Rmed[i-1])*kcoeff[l])/(Rmed[i]-Rmed[i-1]);
      KPhisup = 0.5*(kcoeff[l]+kcoeff[ljp]);
      KPhiinf = 0.5*(kcoeff[l]+kcoeff[ljm]);
      if ( (KRsup < 0) || (KRinf < 0) || (KPhisup < 0) || (KPhiinf < 0) ) {
        printf ("Careful, one of the radiation diffusion coefficients is negative!\n");
        printf ("i = %d: kcoeff_l = %lg, kcoeff_lip = %lg, kcoeff_lim = %lg\n",i,kcoeff[l],kcoeff[lip],kcoeff[lim]);
      }
      /* divF is the divergence of the radiative flux */
      divF = InvDiffRsup[i]*InvRmed[i]*(
					Rsup[i]*KRsup*InvDiffRmed[i+1]*(tempint[lip]-tempint[l]) - \
					Rinf[i]*KRinf*InvDiffRmed[i]*(tempint[l]-tempint[lim]) 
					) +				\
	InvRmed[i]*InvRmed[i]*invdphi*invdphi*(KPhisup*(tempint[ljp]-tempint[l]) - KPhiinf*(tempint[l]-tempint[ljm]));
      raddiff[l] = 2.0*sqrt(ADIABATICINDEX*tempint[l])*pow(Rmed[i],1.5)*divF;
    }
  }
}


int SORsolver_RadiativeDiffusion (Rho, Energy, dt, omega)
     PolarGrid *Rho, *Energy;
     real dt, omega;
{
  /* NEW (Nov 2017): implicit solver for radiative diffusion using
     Successive Over-Relaxation method (adapted from Numerical Recipes
     in C, chap. 19).  we solve dT/dt = div(chi grad(T)) with chi = K
     x 2HK(gamma-1)/sigma = 64sigma_SB lambda
     gamma*(gamma-1)*T^4/Sigma^2/Omega^2/opacity */
  int i, j, l, nr, ns;
  int lip, lim, ljp, ljm;
  real *energy, *tempint, *kcoeff, *chicoeff, *dens, *opacity;
  real *a, *b, *c, *d, *e, *f;
  real bufmsg;
  real dphi, invdphi, lambda, Rl, gradT;
  real chiRsup, chiRinf, chiPhisup, chiPhiinf;
  int n, jchess, ipass;
  real norm_residual,initial_norm_residual=0.0,resid;
  static boolean first=YES;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* Compute Rossland opacities first, used to compute the radiative
     diffusion coefficient 'chi' next */
  ComputeOpacities (Rho, Energy);
  opacity = Opacity->Field;
  tempint = TempInt->Field;
  chicoeff = RadiativeChiCoeff->Field;
  kcoeff = RadiativeKCoeff->Field;
  energy = Energy->Field;
  dens = Rho->Field;
  // Arrays for SOR solver
  a = a_SORarray->Field;
  b = b_SORarray->Field;
  c = c_SORarray->Field;
  d = d_SORarray->Field;
  e = e_SORarray->Field;
  f = f_SORarray->Field;
  //
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  /* First calculate gas temperature */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      tempint[l] = (ADIABATICINDEX-1.0)*energy[l]/dens[l];
    }
  }
  /* Then get chi diffusion coefficient with updated temperature tempint array */
  for (i = 0; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      /* flux-limiting coefficient lambda, as in Kley (1989) */
      gradT = sqrt(pow(InvDiffRmed[i+1]*(tempint[lip]-tempint[l]),2.0) + pow(InvRmed[i]*invdphi*(tempint[ljp]-tempint[l]),2.0));
      Rl = (8.0*sqrt(ADIABATICINDEX*tempint[l])*pow(Rmed[i],1.5)/dens[l]/opacity[l])*(gradT/tempint[l]);
      if (Rl <= 2.0) {
        lambda = 2./(3.+sqrt(9.0+10.0*Rl*Rl));
      } else {
        lambda = 10./(10.*Rl+9.0+sqrt(180.*Rl+81.));
      }
      /* radiative diffusion 'chi' coefficient at i,j (cell-centred) */
      chicoeff[l] = 64.0*sigma_SB*lambda*ADIABATICINDEX*(ADIABATICINDEX-1.0)*pow(tempint[l],4.0)*pow(Rmed[i],3.0)*pow(dens[l],-2.0)/opacity[l];
    }
  }
  // Special case i=nr-1 
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    chicoeff[l] = chicoeff[l-ns];  // I could otherwise do a power-law extrapolation!
  }
  /* Compute array coefficients a to f for SOL solver */
  for (i = 1; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      chiRsup = ((Rmed[i+1]-Rsup[i])*chicoeff[l] + (Rsup[i]-Rmed[i])*chicoeff[lip])/(Rmed[i+1]-Rmed[i]);  // chi_{i+1/2,j}
      chiRinf = ((Rmed[i]-Rinf[i])*chicoeff[lim] + (Rinf[i]-Rmed[i-1])*chicoeff[l])/(Rmed[i]-Rmed[i-1]);  // chi_{i-1/2,j}
      chiPhisup = 0.5*(chicoeff[l]+chicoeff[ljp]); // chi_{i,j+1/2}
      chiPhiinf = 0.5*(chicoeff[l]+chicoeff[ljm]); // chi_{i,j-1/2}
      if ( (chiRsup < 0) || (chiRinf < 0) || (chiPhisup < 0) || (chiPhiinf < 0) ) {
        printf ("Careful, one of the radiation diffusion coefficients is negative in SOL solver!\n");
        printf ("i = %d: chicoeff_l = %lg, chicoeff_lip = %lg, chicoeff_lim = %lg\n",i,chicoeff[l],chicoeff[lip],chicoeff[lim]);
      }
      a[l] = -dt*Rsup[i]*InvDiffRsup[i]*chiRsup*InvRmed[i]*InvDiffRmed[i+1];
      b[l] = -dt*Rinf[i]*InvDiffRsup[i]*chiRinf*InvRmed[i]*InvDiffRmed[i];
      c[l] = -dt*chiPhisup*InvRmed[i]*InvRmed[i]*invdphi*invdphi;
      d[l] = -dt*chiPhiinf*InvRmed[i]*InvRmed[i]*invdphi*invdphi;
      e[l] = 1.0-a[l]-b[l]-c[l]-d[l];
      f[l] = tempint[l];
    }
  }
  /* SOR Solver with Chebyshev acceleration (adapted from Numerical
     Recipes in C, page 869). Greatly inspired from O. Chrenko own
     implementation of SOR solver in his version of Fargo */
  if (first) {
    first = NO;
    if (CPU_Number > 1) {
      ChessBoardIndexing ();
    } else {
      jchess1st = 0;
      jchess2nd = 1;
    }
  }
  // (i) Compute initial norm of residual over active grid of each CPU
  for (i = One_or_active; i < MaxMO_or_active; i++) {
    for (j=0; j<ns; j++) {
      l   = i*ns+j;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      resid = a[l]*tempint[lip] + b[l]*tempint[lim] + c[l]*tempint[ljp] + d[l]*tempint[ljm] + e[l]*tempint[l] - f[l];
      initial_norm_residual += fabs(resid);  // initial norm of residual
    }
  }
  // (ii) Compute initial norm of residual over active grid across all CPUs
  MPI_Allreduce (&initial_norm_residual, &bufmsg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  initial_norm_residual = bufmsg;
  // (iii) Start iterating
  for (n=1; n<=MAXITS; n++) { 
    norm_residual=0.0;  
    for (ipass = 1; ipass <= 2; ipass++) {      // array is swept through as in a chessboard
      if (ipass == 1) jchess = jchess1st;
      if (ipass == 2) jchess = jchess2nd;
      for (i = One_or_active; i < MaxMO_or_active; i++) {
        for (j=jchess; j<ns; j+=2) {    // increment of 'j' index is 2
          l   = i*ns+j;
          lip = l+ns;
          lim = l-ns;
          ljp = l+1;
          if (j == ns-1) ljp = i*ns;
          ljm = l-1;
          if (j == 0) ljm = i*ns+ns-1;
          resid = a[l]*tempint[lip] + b[l]*tempint[lim] + c[l]*tempint[ljp] + d[l]*tempint[ljm] + e[l]*tempint[l] - f[l];
          norm_residual += fabs(resid);
          tempint[l] -= omega*resid/e[l];
        }
        jchess = 1 - jchess; // change the starting value of 'j' for the next increment of 'i'
            }
            if (CPU_Number > 1) SynchronizeOverlapFields (tempint, nr, 1);     // must update the 1st ring of each ghost zone neighbouring to the active mesh
          }
          MPI_Allreduce (&norm_residual, &bufmsg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          norm_residual = bufmsg;
          //masterprint ("n = %d, norm_residual = %lg\n",n,norm_residual);
          // Convergence criterion satisfied:
          if (norm_residual < epsilonSOR*initial_norm_residual) {
            // Now we need to communicate the temperature in the whole set of overlapping rings
            // is that really needed??
            if (CPU_Number > 1) SynchronizeOverlapFields (tempint, nr, CPUOVERLAP);
            // Finally update thermal energy density
            for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          energy[l] = dens[l]*tempint[l] / (ADIABATICINDEX-1.0);
        }
      }
      //masterprint("converged at n=%d as norm_residual = %lg, initial_norm_residual = %lg, ratio=%lg\n",n,norm_residual,initial_norm_residual,norm_residual/initial_norm_residual);
      return n; 
    }
  }
  masterprint("MAXITS exceeded as norm_residual = %lg, initial_norm_residual = %lg, ratio=%lg\n",norm_residual,initial_norm_residual,norm_residual/initial_norm_residual);
  return n;
}

/* Function ensures the odd-even ordering of the SOR method 
 * when the grid is split on multiple CPUs. */
void ChessBoardIndexing ()
/* Function written by O. Chrenko */
{
  int send[3], recv[3], nractive;
  if (CPU_Master) {
    send[0] = 0;        // inner CPU starts with odd cells during the 1st sweep
    send[1] = 1;        // and with even cells during the 2nd sweep
    nractive = NRAD - CPUOVERLAP - 1;   // master has only one ghost zone, but also the innermost ring is skipped in the SOR method
    if (nractive % 2 == 0) {
      send[2] = 0;      // for even number of active rings, next CPU keeps the same indices
    } else {
      send[2] = 1;      // for odd number of active rings, next CPU swaps the indices
    }
    MPI_Send (&send, 3, MPI_INT, CPU_Next, 0, MPI_COMM_WORLD);
  }
  if (CPU_Rank > 0) {
    MPI_Recv (&recv, 3, MPI_INT, CPU_Prev, 0, MPI_COMM_WORLD, &fargostat);
    if (recv[2] == 0) {
      send[0] = recv[0];
      send[1] = recv[1];
    } else {
      send[0] = recv[1];
      send[1] = recv[0];
    }
    if (CPU_Rank != CPU_Highest) {    // outer CPU does not have to send anything
      nractive = NRAD - 2*CPUOVERLAP;   // middle CPUs have two ghost zones
      if (nractive % 2 == 0) {
	send[2] = 0;
      } else {
	send[2] = 1;
      }
      MPI_Send (&send, 3, MPI_INT, CPU_Next, 0, MPI_COMM_WORLD);
    }
  }
  jchess1st = send[0];
  jchess2nd = send[1];
}


int GLOBAL_SORsolver_RadiativeDiffusion (Rho, Energy, dt, omega)
     PolarGrid *Rho, *Energy;
     real dt, omega;
{
  /* NEW (Nov 2017): implicit solver for radiative diffusion using
     Successive Over-Relaxation method (adapted from Numerical Recipes
     in C, chap. 19).  we solve dT/dt = div(chi grad(T)) with chi = K
     x 2HK(gamma-1)/sigma = 64sigma_SB lambda
     gamma*(gamma-1)*T^4/Sigma^2/Omega^2/opacity */
  int i, j, l, nr, ns;
  int lip, lim, ljp, ljm;
  real *energy, *tempint, *kcoeff, *chicoeff, *dens, *opacity;
  real *a, *b, *c, *d, *e, *f;
  real *ga, *gb, *gc, *gd, *ge, *gf, *gtempint;
  real *global_temp;
  real dphi, invdphi, lambda, Rl, gradT;
  real chiRsup, chiRinf, chiPhisup, chiPhiinf;
  int n;
  real norm_residual,initial_norm_residual=0.0,resid;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* Compute Rossland opacities first, used to compute the radiative
     diffusion coefficient 'chi' next */
  ComputeOpacities (Rho, Energy);
  opacity = Opacity->Field;
  tempint = TempInt->Field;
  chicoeff = RadiativeChiCoeff->Field;
  kcoeff = RadiativeKCoeff->Field;
  energy = Energy->Field;
  dens = Rho->Field;
  // Local arrays for SOR solver
  a = a_SORarray->Field;
  b = b_SORarray->Field;
  c = c_SORarray->Field;
  d = d_SORarray->Field;
  e = e_SORarray->Field;
  f = f_SORarray->Field;
  // Global arrays for SOR solver
  ga = Global_a_SORarray->Field;
  gb = Global_b_SORarray->Field;
  gc = Global_c_SORarray->Field;
  gd = Global_d_SORarray->Field;
  ge = Global_e_SORarray->Field;
  gf = Global_f_SORarray->Field;
  gtempint = Global_tempint->Field;
  //
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  /* First calculate gas temperature */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      tempint[l] = (ADIABATICINDEX-1.0)*energy[l]/dens[l];
    }
  }
  /* Then get chi diffusion coefficient with updated temperature tempint array */
  for (i = 0; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      /* flux-limiting coefficient lambda, as in Kley (1989) */
      gradT = sqrt(pow(InvDiffRmed[i+1]*(tempint[lip]-tempint[l]),2.0) + pow(InvRmed[i]*invdphi*(tempint[ljp]-tempint[l]),2.0));
      Rl = (8.0*sqrt(ADIABATICINDEX*tempint[l])*pow(Rmed[i],1.5)/dens[l]/opacity[l])*(gradT/tempint[l]);
      if (Rl <= 2.0) {
        lambda = 2./(3.+sqrt(9.0+10.0*Rl*Rl));
      } else {
        lambda = 10./(10.*Rl+9.0+sqrt(180.*Rl+81.));
      }
      /* radiative diffusion 'chi' coefficient at i,j (cell-centred) */
      chicoeff[l] = 64.0*sigma_SB*lambda*ADIABATICINDEX*(ADIABATICINDEX-1.0)*pow(tempint[l],4.0)*pow(Rmed[i],3.0)*pow(dens[l],-2.0)/opacity[l];
    }
  }
  // Special case i=nr-1 
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    chicoeff[l] = chicoeff[l-ns];  // I could otherwise do a power-law extrapolation!
  }
  /* Compute array coefficients a to f for SOL solver */
  for (i = 1; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      chiRsup = ((Rmed[i+1]-Rsup[i])*chicoeff[l] + (Rsup[i]-Rmed[i])*chicoeff[lip])/(Rmed[i+1]-Rmed[i]);  // chi_{i+1/2,j}
      chiRinf = ((Rmed[i]-Rinf[i])*chicoeff[lim] + (Rinf[i]-Rmed[i-1])*chicoeff[l])/(Rmed[i]-Rmed[i-1]);  // chi_{i-1/2,j}
      chiPhisup = 0.5*(chicoeff[l]+chicoeff[ljp]); // chi_{i,j+1/2}
      chiPhiinf = 0.5*(chicoeff[l]+chicoeff[ljm]); // chi_{i,j-1/2}
      if ( (chiRsup < 0) || (chiRinf < 0) || (chiPhisup < 0) || (chiPhiinf < 0) ) {
	printf ("Careful, one of the radiation diffusion coefficients is negative in SOL solver!\n");
	printf ("i = %d: chicoeff_l = %lg, chicoeff_lip = %lg, chicoeff_lim = %lg\n",i,chicoeff[l],chicoeff[lip],chicoeff[lim]);
      }
      a[l] = -dt*Rsup[i]*InvDiffRsup[i]*chiRsup*InvRmed[i]*InvDiffRmed[i+1];
      b[l] = -dt*Rinf[i]*InvDiffRsup[i]*chiRinf*InvRmed[i]*InvDiffRmed[i];
      c[l] = -dt*chiPhisup*InvRmed[i]*InvRmed[i]*invdphi*invdphi;
      d[l] = -dt*chiPhiinf*InvRmed[i]*InvRmed[i]*invdphi*invdphi;
      e[l] = 1.0-a[l]-b[l]-c[l]-d[l];
      f[l] = tempint[l];
    }
  }
  /* Here we need to define global arrays for tempint, a, b, c, d, e, and f on the whole polar grid! */
  mpi_makeglobalfield (a, ga);
  mpi_makeglobalfield (b, gb);
  mpi_makeglobalfield (c, gc);
  mpi_makeglobalfield (d, gd);
  mpi_makeglobalfield (e, ge);
  mpi_makeglobalfield (f, gf);
  mpi_makeglobalfield (tempint, gtempint);

  /* SOR Solver without Chebyshev acceleration (adapted from Numerical
     Recipes in C, page 869) */
  // (i) Compute initial norm of residual over active grid of each CPU
  for (i = 1; i < GLOBALNRAD-1; i++) {
    for (j = 0; j < ns; j++) {
      l   = i*ns+j;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      resid = ga[l]*gtempint[lip] + gb[l]*gtempint[lim] + gc[l]*gtempint[ljp] + gd[l]*gtempint[ljm] + ge[l]*gtempint[l] - gf[l];
      initial_norm_residual += fabs(resid);  // initial norm of residual
    }
  }
  // (ii) Start iterating
  for (n=1; n<=MAXITS; n++) { 
    norm_residual=0.0;  
    for (i = 1; i < GLOBALNRAD-1; i++) {
      for (j=0; j<ns; j++) {
	l   = i*ns+j;
	lip = l+ns;
	lim = l-ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	resid = ga[l]*gtempint[lip] + gb[l]*gtempint[lim] + gc[l]*gtempint[ljp] + gd[l]*gtempint[ljm] + ge[l]*gtempint[l] - gf[l];
	norm_residual += fabs(resid);
	gtempint[l] -= omega*resid/ge[l];
      }
    }
    //masterprint ("n = %d, norm_residual = %lg\n",n,norm_residual);
    // Convergence criterion satisfied:
    if (norm_residual < epsilonSOR*initial_norm_residual) {
      // Finally update thermal energy density
      for (i = 0; i < nr; i++) {
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  energy[l] = dens[l]*gtempint[(i+IMIN)*NSEC+j] / (ADIABATICINDEX-1.0);
	}
      }
      //masterprint("converged at n=%d as norm_residual = %lg, initial_norm_residual = %lg, ratio=%lg\n",n,norm_residual,initial_norm_residual,norm_residual/initial_norm_residual);
      return n; 
    } 
  }
  masterprint("MAXITS exceeded as norm_residual = %lg, initial_norm_residual = %lg, ratio=%lg\n",norm_residual,initial_norm_residual,norm_residual/initial_norm_residual);
  return n;
}


void SynchronizeOverlapFields (field, nr, nsync)
/* Function originally written by O. Chrenko */
     real *field;
     int nr, nsync;
{
  MPI_Request req1, req2, req3, req4;
  static real *SendBufferInner, *SendBufferOuter;
  static real *RecvBufferInner, *RecvBufferOuter;
  static boolean allocate=YES;
  int prevcpu, nextcpu, size;
  int sendoffsetin, sendoffsetout, recvoffsetin, recvoffsetout;
  /* ----- */
  if (nsync > CPUOVERLAP) {
    printf ("Error! Requested number of rings to be synchronized by the MPI communication is larger than CPUOVERLAP.\n");
    printf ("Terminating now...\n");
    prs_exit (1);
  }
  if (allocate) {
    allocate=NO;
    SendBufferInner = malloc (CPUOVERLAP*NSEC*sizeof(real));    // buffers allocated with the max size
    SendBufferOuter = malloc (CPUOVERLAP*NSEC*sizeof(real));
    RecvBufferInner = malloc (CPUOVERLAP*NSEC*sizeof(real));
    RecvBufferOuter = malloc (CPUOVERLAP*NSEC*sizeof(real));
  }
  size = nsync*NSEC;
  sendoffsetin = CPUOVERLAP*NSEC;
  sendoffsetout = (nr - CPUOVERLAP - nsync)*NSEC;
  recvoffsetin = (CPUOVERLAP - nsync)*NSEC;
  recvoffsetout = (nr - CPUOVERLAP)*NSEC;
  if (CPU_Rank > 0) memcpy (SendBufferInner, field+sendoffsetin, size*sizeof(real));
  if (CPU_Rank != CPU_Highest) memcpy (SendBufferOuter, field+sendoffsetout, size*sizeof(real));
  if (CPU_Rank % 2 == 0) {
    if (CPU_Rank > 0) {
      MPI_Isend (SendBufferInner, size, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Irecv (RecvBufferInner, size, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
    }
    if (CPU_Rank != CPU_Highest) {
      MPI_Isend (SendBufferOuter, size, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
      MPI_Irecv (RecvBufferOuter, size, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
    }
  } else {
    if (CPU_Rank != CPU_Highest) {
      MPI_Irecv (RecvBufferOuter, size, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
      MPI_Isend (SendBufferOuter, size, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
    }
    if (CPU_Rank > 0) {
      MPI_Irecv (RecvBufferInner, size, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Isend (SendBufferInner, size, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
    }
  }
  if (CPU_Rank > 0) {
    MPI_Wait (&req1, &fargostat);
    MPI_Wait (&req2, &fargostat);
    memcpy (field+recvoffsetin, RecvBufferInner, size*sizeof(real));
  }
  if (CPU_Rank != CPU_Highest) {
    MPI_Wait (&req3, &fargostat);
    MPI_Wait (&req4, &fargostat);
    memcpy (field+recvoffsetout, RecvBufferOuter, size*sizeof(real));
  }
}



void ComputeSoundSpeed (Rho, Energy, sys)
     PolarGrid *Rho;
     PolarGrid *Energy;
     PlanetarySystem *sys;
{
  int i, j, l, nr, ns;
  real *dens, *energ, *cs;
  real cs_buf;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    if (!EnergyEquation)
      cs_buf = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!EnergyEquation)
	      cs[l] = cs_buf; // since c_s only depends on radius there
      else
	      cs[l] = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]/dens[l] );
    }
  }
}

/* The pressure and sound speed of the dust's "low-pressure fluid"
shouldn't be given a physical meaning, they are used to avoid
numerical instabilities when handling a fluid with no pressure. That's
why we do not distinguish with or without energy equation... */
void ComputeDustSoundSpeedAndPressure (DRho)
     PolarGrid *DRho;
{
  int i, j, k, l, nr, ns;
  real *drho, *dcs, *dpres, dscbuf;
  nr = DRho->Nrad;
  ns = DRho->Nsec;
  dcs = DSoundSpeed->Field;
  drho = DRho->Field;
  dpres = DPressure->Field;
  for ( i = 0; i < nr; i++ ) {
    dscbuf = DAspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i],DFLARINGINDEX);
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      dcs[l] = dscbuf;
      dpres[l] = drho[l]*dcs[l]*dcs[l];
    }
  }
}

void ComputeDivergenceVelocity (RadialVelocity, AzimuthalVelocity, DRadialVelocity, DAzimuthalVelocity)
     PolarGrid *RadialVelocity, *AzimuthalVelocity;
     PolarGrid *DRadialVelocity, *DAzimuthalVelocity;
{
  int i, j, l, nr, ns;
  int lip, ljp;
  real *vr, *vt, *Dvr, *Dvt, dphi, invdphi;
  real *divergence, *Ddivergence;
  nr  = RadialVelocity->Nrad;
  ns  = RadialVelocity->Nsec;
  vr  = RadialVelocity->Field;
  vt  = AzimuthalVelocity->Field;
  divergence = DivergenceVelocity->Field;
  if (DustFluid) {
    Dvr  = DRadialVelocity->Field;
    Dvt  = DAzimuthalVelocity->Field;
    Ddivergence = DDivergenceVelocity->Field;
  }
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
#pragma omp parallel private(l,lip,ljp,j,ljm,lim)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lip = l+ns;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        divergence[l]  = (vr[lip]*Rsup[i]-vr[l]*Rinf[i])*InvDiffRsup[i]*InvRmed[i];
        divergence[l] += (vt[ljp]-vt[l])*invdphi*InvRmed[i];
        if (DustFluid) {
          Ddivergence[l]  = (Dvr[lip]*Rsup[i]-Dvr[l]*Rinf[i])*InvDiffRsup[i]*InvRmed[i];
          Ddivergence[l] += (Dvt[ljp]-Dvt[l])*invdphi*InvRmed[i];
        }
      }
    }
  }
}
  

void ComputePressureField (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!EnergyEquation) {
        pres[l] = dens[l]*cs[l]*cs[l]; /* since SoundSpeed is not updated */
                                            /* from initialization, cs remains */ 
                                            /* axisymmetric */
      }
      else
	      pres[l] = (ADIABATICINDEX-1.0)*energ[l];
    }
  }
}

void ComputeTemperatureField (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *temp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  temp = Temperature->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!EnergyEquation)
	      temp[l] = MU/R* pres[l]/dens[l];
      else
	      temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
    }
  }
}

real CircumPlanetaryMass (Rho, sys)
     PolarGrid *Rho;
     PlanetarySystem *sys;
{
  int i, j, l, ns;
  real xpl, ypl, rpl;
  real dist, mdcplocal, mdcptotal, MyHillRadius;
  real *dens, *abs, *ord;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  xpl = sys->x[0];
  ypl = sys->y[0];
  rpl = sqrt( xpl*xpl + ypl*ypl );
  mdcplocal = 0.0;
  mdcptotal = 0.0;
  MyHillRadius = rpl * pow( sys->mass[0]/3., 1./3. );
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &fargostat);
  for ( i = Zero_or_active; i < Max_or_active; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      dist = sqrt ( (abs[l]-xpl)*(abs[l]-xpl) +		\
		    (ord[l]-ypl)*(ord[l]-ypl) );
      if ( dist < MyHillRadius ) {
	mdcplocal += Surf[i] * dens[l];
      }
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&mdcplocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    mdcptotal = mdcplocal;
  }
  return mdcptotal;
}


void AddFloorDensity (Rho) 
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  real *dens;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if ( dens[l] < floordens ) 
	dens[l] = floordens;
    }
  }
}


void AddFloorEnergy (Energy) 
     PolarGrid *Energy;
{
  int i, j, l, nr, ns;
  real *energy;
  real e_bottom;
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  energy = Energy->Field;
  for ( i = 0; i < nr; i++ ) {
    e_bottom = 1e-3*EnergyMed[i];  //10^{-3}*initial thermal energy
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if ( energy[l] < e_bottom ) 
	energy[l] = e_bottom;
    }
  }
}


/* Dust velocity in the short-friction time approximation original
   expression from Johansen & Klahr (2005), generalised to include
   dust drag on gas */
void SFTAvelocity (Rho, Vrad, Vtheta, DRho, DVrad, DVtheta)
     PolarGrid *Rho, *DRho, *Vrad, *Vtheta,*DVrad,*DVtheta;
{
  int i, j, l, nr, ns, lim, ljm, lip, ljp;
  real *rho, *press, *drho, *vrad, *vtheta, *dvrad, *dvtheta, *St;
  real gradp_over_rho, dxtheta, invdxtheta, omega, ts;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  press = Pressure->Field;
  drho = DRho->Field;
  dvrad = DVrad->Field;
  dvtheta = DVtheta->Field;
  St = Stokes->Field;

  for (i = 0; i < nr; i++) {
    omega = pow(Rmed[i],-1.5);
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    invdxtheta = 1.0/dxtheta;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ts = St[l]/omega;
      if (DustFeedback)
	      ts *= (1.0+drho[l]/rho[l]);
      lim = l-ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      /* -------------------- */
      /* dust radial velocity */
      /* -------------------- */
      if (i >= 1) {
        //gradp = (cs[l]*cs[l]*rho[l]-cs[lim]*cs[lim]*rho[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
        gradp_over_rho = (press[l]-press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
        //dvrad[l] = vrad[l]+1./sqrt(G*1.0/Rmed[i])*Rmed[i]*St[l]*(gradp+1e-15);
        dvrad[l] = vrad[l] + ts*gradp_over_rho;
      }
      /* ----------------------- */
      /* dust azimuthal velocity */
      /* ----------------------- */
      //gradp = (cs[l]*cs[l]*rho[l]-cs[ljm]*cs[ljm]*rho[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
      gradp_over_rho = (press[l]-press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
      //dvtheta[l] = vtheta[l]+1./sqrt(G*1.0/Rmed[i])*Rmed[i]*St[l]*(gradp+1e-15);
      dvtheta[l] = vtheta[l] + ts*gradp_over_rho;
    }
  }
}


void Diffd (DRho, Rho, DVrad, DVtheta,dt)
     PolarGrid *DRho, *Rho, *DVrad, *DVtheta;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *drho, *rho, *dvrad, *dvtheta, *fdiffrp, *fdifftp, *diag3, *St;
  real dtheta, invdtheta;
  real viscosityp=0.0, viscosity=0.0;
  real frp, frm, ftp, ftm, Dturb;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  
  rho = Rho->Field;
  drho = DRho->Field;
  dvrad = DVrad->Field;
  dvtheta = DVtheta->Field;
  fdiffrp = Fdiffrp->Field;
  fdifftp = Fdifftp->Field;
  diag3 = Diag3->Field;
  St = Stokes->Field;
  // CB: changed the Schmidt number from 1 to (1+4St^2)/[(1+St^2)^2] (St: local Stokes number)

  for (i = 0; i < nr-1; i++){
    dtheta = 2.0*PI/(real)ns;
    invdtheta = 1.0/dtheta;
    if ((ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0)) {
      viscosityp = FViscosity (Rsup[i]);
      viscosity = FViscosity (Rmed[i]);
    }
    for (j = 0; j < ns; j++){
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      Dturb = viscosityp*(1.0+4.0*St[lip]*St[lip])*pow((1.0+St[lip]*St[lip]),-2.0);
      fdiffrp[l]=-Dturb*(rho[lip]+rho[l])/2.*(drho[lip]/rho[lip]-drho[l]/rho[l])*InvDiffRmed[i+1];
      Dturb = viscosity*(1.0+4.0*St[l]*St[l])*pow((1.0+St[l]*St[l]),-2.0);
      fdifftp[l]=-Dturb*(rho[ljp]+rho[l])/2.*(drho[ljp]/rho[ljp]-drho[l]/rho[l])*invdtheta/Rmed[i];
      Dturb = viscosityp*(1.0+4.0*St[lip]*St[lip])*pow((1.0+St[lip]*St[lip]),-2.0);
      diag3[l]=-Dturb*(log(drho[lip]/rho[lip])-log(drho[l]/rho[l]))/(log(Rmed[i+1])-log(Rmed[i]))/Rsup[i];
    }
  } 
  for (i = 1; i < nr-1; i++){
    if(Rmed[i]>DUSTDIFFINN){
      for (j = 0; j < ns; j++){
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        frp=fdiffrp[l];
        frm=fdiffrp[lim];
        ftp=fdifftp[l];
        ftm=fdifftp[ljm];
        drho[l]=drho[l]-(Rsup[i]*dtheta*frp-Rinf[i]*dtheta*frm+(Rsup[i]-Rinf[i])*(ftp-ftm))*dt*InvSurf[i];
      }
    }
  }
}
