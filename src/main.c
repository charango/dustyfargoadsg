/** \file main.c

Main file of the distribution. Manages the call to initialization
functions, then the main loop.

*/

#include "mp.h"

boolean         TimeToWrite, Restart = NO, OpenInner = NO, OpenInnerDust = NO;
int             begin_i = 0, NbRestart = 0, verbose = NO;
static int      InnerOutputCounter=0, StillWriteOneOutput;
extern real     LostMass,LostMassd,AccRate,AccRated;
extern boolean  Corotating, ReadPlanetFileAtRestart, Write_DustDensity, Evanescent, DampToViscous;
extern boolean  CorotateWithOuterPlanet, DiscEvaporation, Write_Jacobi, RestartWithNewDust, ComputeCPDMass;
extern boolean  SelfGravity, SGZeroMode, EnergyEquation, SoftWriting, Write_DustSystem, DustFluid, DustFeedback;
extern boolean  CompareSGAndSummationTorques;
real            ScalingFactor = 1.0;

int
main(argc, argv)
     int argc;
     char *argv[];
{
  PolarGrid   *gas_density;
  PolarGrid   *gas_v_rad; 
  PolarGrid   *gas_v_theta; 
  PolarGrid   *gas_energy; 
  PolarGrid   *gas_label;
  PolarGrid   *dust_pc_density;
  PolarGrid   *dust_density;
  PolarGrid   *dust_v_rad;
  PolarGrid   *dust_v_theta;
  int          i,j,k;
  real         foostep = 0.;
  real         r, v;
  real         Initial_Disc_Dust_Mass;
  boolean      disable = NO, TimeInfo = NO, Profiling = NO;
  boolean      Stockholm = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;
  Force *force;
  DustSystem *dustsys;
  real axidens_idm[GLOBALNRAD];
  FILE *fich_idm;
  char name_idm[256];


  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  setfpe ();  /* Control behavior for floating point
		 exceptions trapping (default is not to do anything) */
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-secndovtpfamzib0123456789") != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 'n'))
	disable = YES;
      if (strchr (argv[i], 'e'))
	Stockholm = YES;
      if (strchr (argv[i], 'v'))
	verbose = YES;
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'c'))
	SloppyCFL = YES;
      if (strchr (argv[i], 'p'))
	Profiling = YES;
      if (strchr (argv[i], 'd'))
	debug = YES;
      if (strchr (argv[i], 'b'))
	CentrifugalBalance = YES;
      if (strchr (argv[i], 'm'))
	Merge = YES;
      if (strchr (argv[i], 'a'))
	MonitorIntegral = YES;
      if (strchr (argv[i], 'z'))
	FakeSequential = YES;
      if (strchr (argv[i], 'i')) {
	StoreSigma = YES;
	if (EnergyEquation)
	  StoreEnergy = YES;
      }
      if (strchr (argv[i], '0'))
	OnlyInit = YES;
      if ((argv[i][1] >= '1') && (argv[i][1] <= '9')) {
	GotoNextOutput = YES;
	StillWriteOneOutput = (int)(argv[i][1]-'0');
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'o')) {
	OverridesOutputdir = YES;
	i++;
	sprintf (NewOutputdir, "%s", argv[i]);
      } else {
	if (strchr (argv[i], 'f')) {
	  i++;
	  ScalingFactor = atof(argv[i]);
	  masterprint ("Scaling factor = %g\n", ScalingFactor);
	  if ((ScalingFactor <= 0)) {
	    masterprint ("Incorrect scaling factor\n");
	    PrintUsage (argv[0]);
	  }
	}
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }
  
  if ( (StoreSigma || StoreEnergy) && !(Restart)) {
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("or surface internal energy in a non-restart run.\n");
    mastererr ("Aborted\n");
    prs_exit (0);
  }
  
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);
  ReadVariables (ParameterFile);
  SplitDomain ();
  
  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    prs_exit (0);
  
  DumpSources (argc, argv, ParameterFile, PLANETCONFIG);
  ComputeCodeUnits ();
  
  masterprint ("Allocating arrays...");
  fflush (stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "gasdens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "gasvrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "gasvtheta");
  gas_energy         = CreatePolarGrid(NRAD, NSEC, "gasenergy");
  gas_label          = CreatePolarGrid(NRAD, NSEC, "gaslabel");
  dust_pc_density    = CreatePolarGrid(NRAD, NSEC, "pcdens");
  dust_density       = CreatePolarGrid(NRAD, NSEC, "dustdens");
  dust_v_rad         = CreatePolarGrid(NRAD, NSEC, "dustvrad");
  dust_v_theta       = CreatePolarGrid(NRAD, NSEC, "dustvtheta");
  torquesumdisc      = CreatePolarGrid(NRAD, NSEC, "torquesumdisc");
  masterprint ("done.\n");
  FillPolar1DArrays ();
  force = AllocateForce ();
  
  /* Here planets are initialized feeling star potential but they do
     not feel disk potential */
  sys = InitPlanetarySystem (PLANETCONFIG);
  
  /* Gas density initialization */
  InitGasDensity (gas_density);

  /* Dust density initialization */
  if (DustFluid)
    InitDustDensity (dust_density);

  /* If energy equation is taken into account, we initialize the gas
     thermal energy */
  if (EnergyEquation)
    InitGasEnergy (gas_energy);
  
  if ( SelfGravity ) {
    /* If SelfGravity = YES or Z, planets are initialized feeling disk
       potential. Only the surface density is required to calculate
       the radial self-gravity acceleration. The disk radial and
       azimutal velocities are not updated. If this is a restart run,
       we first need to read the density at restart! */
    if (Restart == YES)
      ReadfromFile (gas_density, "gasdens", NbRestart);
    // compute self-gravity w/o updating gas velocities
    compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, dust_v_rad, dust_v_theta, foostep, NO);
    init_planetarysys_withSG (sys);
  }

  /* Stores the mass of the super-particles for calculating the
     dust-to-gas density ratio */
  if (NBPART != 0) {
    if ( (Restart == NO) || ( (Restart == YES) && (RestartWithNewDust == YES) ) ) {
      // initial dust surface density profile known by all CPUs
      mpi_make1Dprofile (dust_density->Field, GLOBAL_bufarray);
      Initial_Disc_Dust_Mass = 0.0;
      // below we compute the total dust disc mass between RMINDUST and RMAXDUST
      for (i=0; i<GLOBALNRAD; i++) {
	if ( (Radii[i] >= RMINDUST) && (Radii[i] <= RMAXDUST) )
	  Initial_Disc_Dust_Mass += (GLOBAL_bufarray[i]*0.5*(PMAX-PMIN)*(Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]));
      }
      Particles_Mass_Initial = Initial_Disc_Dust_Mass/NBPART;
      Particles_Mass = Particles_Mass_Initial;
      // Keep track of disc's initial mass in a file for a restart run
      if (CPU_Master) {
	sprintf (name_idm, "%s%s.dat", OUTPUTDIR, "particlesmass");
	fich_idm = fopen(name_idm, "w");
	if (fich_idm == NULL) {
	  fprintf (stderr, "Can't write 'particlesmass.dat' file. Aborting.\n");
	  prs_exit (1);
	}
	fprintf (fich_idm, "%lg",Particles_Mass);
	fclose (fich_idm);
      }
    } else {
      if (DustFeedback) {
	// Case of a restart simulation with already evolved dust particles
	sprintf (name_idm, "%s%s.dat", OUTPUTDIR, "particlesmass");
	fich_idm = fopen(name_idm, "r");
	if (fich_idm == NULL) {
	  fprintf (stderr, "Can't read 'particlesmass.dat' file. Aborting.\n");
	  prs_exit (1);
	}
	fscanf (fich_idm, "%lg",&Particles_Mass);
	Particles_Mass_Initial = Particles_Mass;
      }
    }
    masterprint ("Mass of the super-particles is %lg\n",Particles_Mass);
  }

  ListPlanets (sys);
  OmegaFrame = OMEGAFRAME;
  if (Corotating == YES) {
    if (!SelfGravity)
      OmegaFrame = GetPsysInfo (sys, FREQUENCY);
    else {
      if (!CorotateWithOuterPlanet) {
	r = sqrt( sys->x[0]*sys->x[0] + sys->y[0]*sys->y[0] );
	v = sqrt( sys->vx[0]*sys->vx[0] + sys->vy[0]*sys->vy[0] );
      } else {
	r = sqrt( sys->x[1]*sys->x[1] + sys->y[1]*sys->y[1] );
	v = sqrt( sys->vx[1]*sys->vx[1] + sys->vy[1]*sys->vy[1] );
      }
      /* If self-gravity is included, OmegaFrame is first calculated
	 as the angular frequency the planet would have if it was on a
	 fixed circular orbit at x=sys->x[0], with initial velocity
	 sys->vy[0], which includes the radial initial
	 self-gravitating acceleration */
      OmegaFrame = v / r;
    }
  }

  /* Only gas and dust velocities remain to be initialized */
  Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys, dust_density, dust_v_rad, dust_v_theta);
  
  /* as of May 2016: evanescent boundary condition can now damp toward
     1D viscously evolving density and radial velocity profiles. Both
     profiles are solved through an independent 1D viscous evolution
     problem without planets. */
  if ( (Evanescent && DampToViscous) || (DampToViscous))
    InitializeOneDViscousEvolution ();

  /* Initial gas_density is used to compute the circumplanetary mass
     with initial density field */
  if (ComputeCPDMass)
    mdcp0 = CircumPlanetaryMass (gas_density, sys);
    
  if (Restart == YES) {
    begin_i         = NbRestart * NINTERM;
    if (ReadPlanetFileAtRestart)
      RestartPlanetarySystem (NbRestart, sys);
    LostMass = GetfromPlanetFile (NbRestart, 7, 0); /* 0 refers to planet #0 */
    PhysicalTime  = GetfromPlanetFile (NbRestart, 8, 0);
    if (Corotating == YES)
      OmegaFrame = GetfromPlanetFile (NbRestart, 9, 0); // otherwise OmegaFrame will be 0
  } else {			/* We initialize 'planet[i].dat' file */
    EmptyPlanetSystemFile (sys);
  }
  
  /* Here we initialize the dust as Lagrangian particles */
  if (NBPART != 0) {
    dustsys = InitDustSystem ();
    if (Write_Jacobi) {
      for (k=0; k<NBPART; k++)
	CreateJacobiFile(dustsys,k);
    }
    EmptyDustSystemFile (dustsys, TimeStep);    /*We initialize 'Dustsystat[i]' file*/
    /* at this stage, all cpus know about all the dust particles -- 
       call to interpolation helps initialize dust pc density array as 
       well as dust feedback arrays */
    interpolation(dustsys, gas_v_rad, gas_v_theta, gas_density, dust_pc_density, TimeStep);
  }
  
  if (MonitorIntegral == YES)
    CheckMomentumConservation (gas_density, gas_v_theta, sys);
  
  PhysicalTimeInitial = PhysicalTime;
  
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);
  MultiplyPolarGridbyConstant (dust_density, ScalingFactor); /* Polar grid dust scaling, used with -f options during start/restart */

  
  /* ========================== */
  /* LOOP OVER TIME STARTS HERE */
  /* ========================== */
  for (i = begin_i; i <= NTOT; i++) {
    InnerOutputCounter++;
    
    /* SoftWriting boolean is false by default. Default outputs for
       tqwk*.dat and bigplanet*.dat files */
    if ( (!SoftWriting) && (InnerOutputCounter == 1)) {
      InnerOutputCounter = 0;
      WriteBigPlanetSystemFile (sys, TimeStep);
      SolveOrbits (sys);
      UpdateLog (force, sys, gas_density, TimeStep, PhysicalTime);
      if (Stockholm == YES)
	UpdateLogStockholm (sys, gas_density, TimeStep, PhysicalTime);
      if (SelfGravity) {
	// we just compute self-gravitating acceleration prior to computing SG-interpolated torque on planet
	compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, dust_v_rad, dust_v_theta, foostep, YES);
	if (!SGZeroMode)
	  UpdateLogSG (sys, TimeStep, PhysicalTime);
      }
    }

     /* Outputs every NINTERM timesteps are done here */
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) {
      TimeToWrite = YES;

      /* NEW (June 2023) */
      if (CompareSGAndSummationTorques) {
	// compute self-gravity w/o updating gas velocities
	if (SelfGravity)
	  compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, dust_v_rad, dust_v_theta, foostep, YES);
	CompareSGandSummationTorques (force, gas_density, sys);
      }
      
      SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, dust_pc_density, dust_density, dust_v_rad, dust_v_theta, dustsys);
      WritePlanetSystemFile (sys, TimeStep);
      if (Evanescent && DampToViscous)
	Write1DViscProfiles (TimeStep);
      if (NBPART != 0) {
	if (Write_DustSystem)
	  WriteDustSystemFile (dustsys, TimeStep);
	if (Write_Jacobi) {
	  for(j=0; j<NBPART; j++)
	    WriteJacobi(dustsys,sys,j);
	}
      }
      
      /* If SoftWriting true, tqwk* files and bigplanet*.dat files are
	 written every NINTERM timesteps */
      if (SoftWriting) {
	WriteBigPlanetSystemFile (sys, TimeStep);
	SolveOrbits (sys);
	UpdateLog (force, sys, gas_density, TimeStep, PhysicalTime);
	if (Stockholm == YES)
	  UpdateLogStockholm (sys, gas_density, TimeStep, PhysicalTime);
	if (SelfGravity) {
	  // we just compute self-gravitating acceleration prior to computing SG-interpolated torque on planet
	  compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, dust_v_rad, dust_v_theta, foostep, YES);
	  if (!SGZeroMode)
	    UpdateLogSG (sys, TimeStep, PhysicalTime);
	}
      }
      
      if ((OnlyInit) || ((GotoNextOutput) && (!StillWriteOneOutput))) {
	MPI_Finalize();
	return 0;
      }
      
      StillWriteOneOutput--;
      if (TimeInfo == YES)	/* Time monitoring is done here */
	GiveTimeInfo (TimeStep);
      
    }
    else {
      TimeToWrite = NO;
    }
    
    /***********************/
    /* Hydrodynamical Part */
    /***********************/
    InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
    AlgoGas (force, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, dust_density, dust_pc_density, dust_v_rad, dust_v_theta, sys, dustsys);
    GiveSpecificTime (Profiling, t_Hydro);
    
    if (MonitorIntegral == YES) {
      CheckMomentumConservation (gas_density, gas_v_theta, sys);
      masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
      masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
      masterprint ("Gas total Energy : %.18g\n", GasTotalEnergy (gas_density, gas_v_rad, gas_v_theta, gas_energy));
    }
    
  } // end loop for (i = begin_i; i <= NTOT; i++)
  
  if (CPU_Master)
    printf("Dt min = %.12g and Dt max = %.12g at the end of the simulation\n", MINDT, MAXDT);
  
  FreePlanetary (sys);
  FreeForce (force);
  if (NBPART != 0)
    FreeDust (dustsys);
  if ( SelfGravity && !SGZeroMode ) {
    rfftwnd_mpi_destroy_plan(SGP_fftplan_forward);
    rfftwnd_mpi_destroy_plan(SGP_fftplan_backward);
  }
  
  MPI_Finalize ();
  return 0;
}
