/** \file merge.c

Contains the function that merges the output of different processors.
The resulting merged file is undistinguishable from the file that
would have been produced by a sequential run.

*/

#include "mp.h"

void merge (nb)
     int nb;
{
  extern boolean  Write_Density, Write_Velocity, Write_Energy;
  extern boolean  Write_Temperature, Write_DivV, Write_ViscHeat, Write_TherDiff, Write_RadDiff, Write_TherCool;
  extern boolean  Write_Potential, Write_Test, Write_gr, Write_gtheta, Write_DustDensity, Write_RadFBAcc, Write_AziFBAcc; 
  extern boolean SelfGravity, SGZeroMode;
  extern boolean  AdvecteLabel, DustFluid;
  boolean bool = NO;
  int i, j, l;
  int one_if_odd;
  char radix[512];
  char command[1024];
  if (!CPU_Master) return;
  message ("Merging output files...");
  for (j = 0; j < 19+(AdvecteLabel == YES); j++) {
    switch (j) {
    case 0:
      strcpy (radix, "gasdens");
      bool = Write_Density;
      break;
    case 1: 
      strcpy (radix, "gasvrad");
      bool = Write_Velocity;
      break;
    case 2: 
      strcpy (radix, "gasvtheta");
      bool = Write_Velocity;
      break;
    case 3: 
      strcpy (radix, "gasenergy");
      bool = Write_Energy;
      break;
    case 4: 
      strcpy (radix, "Temperature");
      bool = Write_Temperature;
      break;
    case 5: 
      strcpy (radix, "DivV");
      bool = Write_DivV;
      break;
    case 6: 
      strcpy (radix, "ViscousHeating");
      bool = Write_ViscHeat;
      break;
    case 7: 
      strcpy (radix, "EntropyDiff");
      bool = Write_TherDiff;
      break;
    case 8: 
      strcpy (radix, "RadiativeKCoeff");
      bool = Write_RadDiff;
      break;
    case 9: 
      strcpy (radix, "ThermalCooling");
      bool = Write_TherCool;
      break;
    case 10: 
      strcpy (radix, "Potential");
      bool = Write_Potential;
      break;
    case 11: 
      strcpy (radix, "Test");
      bool = Write_Test;
      break;
    case 12: 
      strcpy (radix, "sgaccr");
      bool = Write_gr;
      break;
    case 13: 
      strcpy (radix, "sgacctheta");
      bool = Write_gtheta;
      break;
    case 14: 
      strcpy (radix, "torquesg");
      bool = Write_gtheta;
      break;
    case 15: 
      strcpy (radix, "torquesumdisc");
      bool = Write_gtheta;
      break;
    case 16: 
      strcpy (radix, "dust");
      bool = Write_DustDensity*(NBPART == 0 ? 0 : 1);
      break;
    case 17: 
      strcpy (radix, "RadFBAcc");
      bool = Write_RadFBAcc;
      break;
    case 18: 
      strcpy (radix, "AziFBAcc");
      bool = Write_AziFBAcc;
      break;
    case 19: 
      strcpy (radix, "label");
      bool = YES;
      break;
    }
    one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
    if (bool == YES) {
      if ( SelfGravity && !SGZeroMode ) {
	for (i = 0; i < (CPU_Number+one_if_odd)/2; i++) {
	  if ( i != 0 ) {
	    sprintf (command, "cd %s; cat %s%d.dat.%05d >> %s%d.dat", \
		     OUTPUTDIR, radix, nb, i, radix, nb);
	    system (command);
	  }
	  l = (i + (CPU_Number + one_if_odd)/2)%(CPU_Number + one_if_odd);
	  if ( i != CPU_Highest ) {
	    sprintf (command, "cd %s; cat %s%d.dat.%05d >> %s%d.dat", \
		     OUTPUTDIR, radix, nb, l, radix, nb);
	    system (command);
	  }
	}
	sprintf (command, "cd %s; rm -f %s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
      else {
	for (i = 1; i < CPU_Number; i++) {
	  sprintf (command, "cd %s; cat %s%d.dat.%05d >> %s%d.dat", \
		   OUTPUTDIR, radix, nb, i, radix, nb);
	  system (command);
	}
	sprintf (command, "cd %s; rm -f %s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
    }
  }
  message ("done\n");
  fflush (stdout);
}


void mergedust (nb)
     int nb;
{
  extern boolean Write_Density, Write_DustDensity, Write_Velocity, Write_StokesNumber;
  extern boolean SelfGravity, SGZeroMode;
  boolean bool = NO;
  int i, j, l;
  int one_if_odd;
  char radix[512];
  char command[1024];
  if (!CPU_Master) return;
  message ("Merging output files...");
  for (j = 0; j < 3 + (Write_StokesNumber == YES ? 1 : 0); j++) {
    switch (j) {
    case 0:
      strcpy (radix, "dustdens");
      bool = Write_Density;
      break;
    case 1: 
      strcpy (radix, "dustvrad");
      bool = Write_Velocity;
      break;
    case 2: 
      strcpy (radix, "dustvtheta");
      bool = Write_Velocity;
      break;
    case 3:
      strcpy (radix, "Stokes");
      bool = Write_StokesNumber;
      break;
    }
    one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
    if (bool == YES) {
      if ( SelfGravity && !SGZeroMode ) {
	for (i = 0; i < (CPU_Number+one_if_odd)/2; i++) {
	  if ( i != 0 ) {
	    sprintf (command, "cd %s; cat %s%d.dat.%05d >> %s%d.dat", \
		     OUTPUTDIR, radix, nb, i, radix, nb);
	    system (command);
	  }
	  l = (i + (CPU_Number + one_if_odd)/2)%(CPU_Number + one_if_odd);
	  if ( i != CPU_Highest ) {
	    sprintf (command, "cd %s; cat %s%d.dat.%05d >> %s%d.dat", \
		     OUTPUTDIR, radix, nb, l, radix, nb);
	    system (command);
	  }
	}
	sprintf (command, "cd %s; rm -f %s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
      else {
	for (i = 1; i < CPU_Number; i++) {
	  sprintf (command, "cd %s; cat %s%d.dat.%05d >> %s%d.dat", \
		   OUTPUTDIR, radix, nb, i, radix, nb);
	  system (command);
	}
	sprintf (command, "cd %s; rm -f %s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
    }
  }
  message ("done\n");
  fflush (stdout);
}
