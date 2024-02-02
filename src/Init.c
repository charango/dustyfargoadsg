/** \file Init.c

Contains the functions needed to initialize the hydrodynamics arrays.
These can be initialized by reading a given output (in the case of a
restart) or by calling a function, InitEuler (), which contains
analytic prescription for the different hydrodynamics fields. Note
that this function InitEuler() is located in SourceEuler.c, which
itself calls InitGas(), in the file Pframeforce.c.
Also, note that the present file contains InitLabel(), which sets
the initial value of a passive scalar.
*/

#include "mp.h"

extern boolean Restart;
extern int     NbRestart;

void ReadfromFile (array, fileprefix, filenumber)
     PolarGrid *array;
     char *fileprefix;
     int filenumber;
{
  int nr, ns, c, foo=0;
  real *field;
  char name[256];
  FILE *input;
  /* Simultaneous read access to the same file have been observed to
     give wrong results. */
  /* A sequential reading is imposed below. */
  /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Prev, 10, MPI_COMM_WORLD, &fargostat);
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
    return;
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); 
    /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message
     that it expects */
  if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything
				   meanwhile */
}

void InitLabel (array, sys)
     PolarGrid *array;
     PlanetarySystem *sys;
{
  int nr, ns, i, j, l;
  real xp, yp, rp;
  real x, y, angle, distance, rhill;
  real *field;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  xp = sys->x[0];
  yp = sys->y[0];
  rp = sqrt ( xp*xp + yp*yp );
  rhill = rp * pow( sys->mass[0]/3., 1./3 );
   /* Initialize label as you wish. In this example, label only takes
      into account fluid elements inside the planet's Hill Sphere */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      if ( (Rmed[i] > 1.5) && (Rmed[i] < 2.3) )
	field[l] = 1.;
      else 
	field[l] = 0.;
      /*
      angle = Azimuth[j];
      x = Rmed[i] * cos(angle);
      y = Rmed[i] * sin(angle);
      distance = sqrt( (x-xp)*(x-xp) + (y-yp)*(y-yp) );
      if ( distance < rhill )
	field[l] = 1.0;
      else
	field[l] = 0.0;
      */
    }
  }
}

void Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, pla_sys, dust_density, dust_v_rad, dust_v_theta)
     PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_energy, *gas_label;
     PolarGrid *dust_density, *dust_v_rad, *dust_v_theta;
     PlanetarySystem *pla_sys;
     /* CB: NEED TO ADD INITIALIZATION OF DUST FLUID */
{
  extern boolean EnergyEquation, EntropyDiffusion, RadiativeDiffusion, ImplicitRadiativeDiffusion, ThermalCooling, DustFluid, RestartWithNewDust;
  real *energ, *dens;
  FILE *output;
  char OutputName[256];
  int i, j, l, nr, ns;
  energ = gas_energy->Field;
  nr = gas_energy->Nrad;
  ns = gas_energy->Nsec;
  ReadPrevDim ();
  InitEuler (gas_v_rad, gas_v_theta, gas_density, gas_energy, dust_v_rad, dust_v_theta, dust_density, pla_sys);
  InitLabel (gas_label, pla_sys);
  if (Restart == YES) {
    CheckRebin (0);
    CheckRebin (NbRestart);
    /* Now that OldRmed is built, Cpu_master writes again Radii in file
     used_rad.dat */
    if ( CPU_Master ) {
      sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
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
    MPI_Barrier (MPI_COMM_WORLD);
    /* Don't start reading before master has finished rebining... */
    /* It shouldn't be a problem though since a sequential read is */
    /* imposed in the ReadfromFile function below */
    mastererr ("Reading restart files...");
    fflush (stderr);
    ReadfromFile (gas_density, "gasdens", NbRestart);
    ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
    ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
    if (EnergyEquation) {
      ReadfromFile (gas_energy, "Temperature", NbRestart);
      /* ! gas_energy accounts for the gas temperature... */
      dens = gas_density->Field;
      for (i=0; i<nr; i++) {
	for (j=0; j<ns; j++) {
	  l = i*ns + j;
	  energ[l] = dens[l]*energ[l]/(ADIABATICINDEX-1.0);
	  /* this is e = dens*temp / (gamma-1) */
	}
      }
    }
    if (DustFluid && (!RestartWithNewDust)) {
      ReadfromFile (dust_density, "dustdens", NbRestart);
      ReadfromFile (dust_v_rad, "dustvrad", NbRestart);
      ReadfromFile (dust_v_theta, "dustvtheta", NbRestart);
    }
    /* To output restart fields correctly */
    ComputeSoundSpeed (gas_density, gas_energy, pla_sys);
    ComputePressureField (gas_density, gas_energy);
    ComputeTemperatureField (gas_density, gas_energy);
    if (EntropyDiffusion)
      ComputeEntropyDiffusion (gas_density, gas_energy);
    if (RadiativeDiffusion && !ImplicitRadiativeDiffusion)
      ComputeRadiativeDiffusion (gas_density, gas_energy);
    if (ThermalCooling)
      ComputeThermalCooling (gas_density, gas_energy);
    ComputeViscousTerms (gas_v_rad, gas_v_theta, gas_density,dust_v_rad, dust_v_theta, dust_density);
    ComputeViscousHeating (gas_density);
    ReadfromFile (gas_label, "gaslabel", NbRestart);
    if (StoreSigma) RefillSigma (gas_density);
    if (StoreSigma && DustFluid) RefillDSigma (dust_density);
    /* StoreEnergy = NO if EnergyEquation = NO */
    if (StoreEnergy) RefillEnergy (gas_energy);
    fprintf (stderr, "done\n");
    fflush (stderr);
  }
  WriteDim (); 
}
