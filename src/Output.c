/** \file Output.c

Contains most of the functions that write the output files.  In
addition to the writing of hydrodynamics files (handled by SendOutput
()), this file also contains the functions that update the planet.dat
and bigplanet.dat files, and the functions that seek information about
the planets at a restart.
*/

#include "mp.h"

static real     Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual;
extern real     LostMass;
extern boolean  Write_Density, Write_Velocity, Write_Energy, IsDisk;
extern boolean  Write_Temperature, Write_DivV, Write_TherDiff, Write_RadDiff, Write_TherCool, Write_ViscHeat;
extern boolean  Write_Potential, Write_Test, Write_DustDensity, Write_StokesNumber, Write_RadFBAcc, Write_AziFBAcc;
extern boolean  Write_gr, Write_gtheta, DustFeedback;
extern boolean  AdvecteLabel;

void EmptyPlanetSystemFile (sys)
     PlanetarySystem *sys;
{
  FILE *output;
  char name[256];
  int i, n;
  n = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < n; i++) {
    sprintf (name, "%splanet%d.dat", OUTPUTDIR, i);
    output = fopen (name, "w");
    if (output == NULL) {
      fprintf (stderr, "Can't write %s file. Aborting.\n", name);
      prs_exit (1);
    }
    fclose (output);
  }
}

void WritePlanetFile (timestep, n)
     int timestep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  printf ("Updating 'planet%d.dat'...", n);
  fflush (stdout);
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'planet%d.dat' file. Aborting.\n", n);
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\n", timestep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WritePlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WritePlanetFile (t, i);
  }
}
   

void WriteBigPlanetFile (timestep, n)
     int timestep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  sprintf (name, "%sbigplanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'bigplanet.dat' file. Aborting.\n");
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\n", timestep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp);
  fclose (output);
}

void WriteBigPlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WriteBigPlanetFile (t, i);
  }
}

real GetfromPlanetFile (timestep, column, n)
     int timestep, column, n;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'planet%d.dat' file. Aborting restart.\n",n);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'planet%d.dat'. Aborting restart.\n",n);
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != timestep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'planet%d.dat' file. Aborting restart.\n", timestep,n);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}

void RestartPlanetarySystem (timestep, sys)
     PlanetarySystem *sys;
     int timestep;
{
  int k;
  real bufmass, bufvx, bufvy;
  for (k = 0; k < sys->nb; k++) {
    sys->x[k] = GetfromPlanetFile (timestep, 2, k);
    sys->y[k] = GetfromPlanetFile (timestep, 3, k);
    bufvx   = GetfromPlanetFile (timestep, 4, k);
    bufvy   = GetfromPlanetFile (timestep, 5, k);
    bufmass = GetfromPlanetFile (timestep, 6, k);
    InitialPlanetMass[k] = bufmass;
    FinalPlanetMass[k] = bufmass;
    /*
      if ( fabs((bufmass-sys->mass[k])/bufmass) > 1e-4 ) {
      masterprint ("There is a discrepancy between the planet mass of your restarting calculation (mass in planet0.dat, at the restarting line), and the planet mass you mention in the .cfg file. I guess you know what you do, so I will consider the planet mass in your .cfg restarting file. The planet velocities have been recalculated with the new planet mass so as to yield a rotational equilibrium. Please check this is really what you meant to do!\n");
      sys->vx[k] = bufvx*sqrt( (1.+sys->mass[k])/(1.+bufmass) );
      sys->vy[k] = bufvy*sqrt( (1.+sys->mass[k])/(1.+bufmass) );
      }
      else {
    */
    sys->vx[k] = bufvx;
    sys->vy[k] = bufvy;
    sys->mass[k] = bufmass;
    masterprint ("At restart, mass of planet %d is %.15g, its y is %lg\n", k, sys->mass[k], sys->y[k]);
    //}
  }
}

void WriteDiskPolar(array, number)   // Ecriture dans un fichier a partir d'un tableau DEJA CONNU !
     PolarGrid 	*array;
     int 	number;
{
  int           Nr, Ns;
  FILE          *dump;
  char 		name[256];
  real 		*ptr;
  ptr = array->Field;
  if (CPU_Master)
    sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Nr = array->Nrad;
  Ns = array->Nsec;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s%d.dat'...", array->Name, number);
  fflush (stdout);
  MPI_Barrier (MPI_COMM_WORLD);
/* We strip the first CPUOVERLAP rings if the current CPU is not the 
   innermost one */
  if (CPU_Rank > 0) {
    ptr += CPUOVERLAP*Ns;
    Nr -=CPUOVERLAP ;
  }
/* We strip the last CPUOVERLAP rings if the current CPU is not the outermost
   one, equal to CPU_Highest in all cases */
  if (CPU_Rank != CPU_Highest) {
    Nr -=CPUOVERLAP;
  }
  fwrite (ptr, sizeof(real), Nr*Ns,dump); //on écrit ptr dans dump.. pb : d'où vient le "array" que l'on a mis dans ptr ??
  fclose(dump);
  //fprintf(stdout, "%d/", CPU_Rank);  
  fflush(stdout);
  MPI_Barrier (MPI_COMM_WORLD);
  masterprint("done\n");
}

void WriteDim () {	  
  char filename[256];
  FILE 	*dim;
  if (!CPU_Master) return;
  sprintf (filename, "%sdims.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    prs_exit (1);
  }
  fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n",0,0,0,0,RMAX, NTOT/NINTERM, GLOBALNRAD, NSEC);
  fclose (dim);
}


void SendOutput (index, dens, gasvr, gasvt, gasenerg, label, dustpcdens, dustdens, dustvr, dustvt, dustsys)
     int          index;
     PolarGrid   *dens, *gasvr, *gasvt, *label, *gasenerg, *dustpcdens;
     PolarGrid   *dustdens, *dustvr, *dustvt;
     DustSystem  *dustsys;
{
  extern boolean DustFluid, SelfGravity;
  if (CPU_Master)
    printf ("\n*** OUTPUT %d ***\n", index);
  if (IsDisk == YES) {
    if (AdvecteLabel == YES) WriteDiskPolar (label, index);
    if (Write_Density == YES) {
      WriteDiskPolar (dens, index);
      if (DustFluid)
	WriteDiskPolar (dustdens, index);
    }
    if (Write_Velocity == YES) {
      WriteDiskPolar (gasvr, index);
      WriteDiskPolar (gasvt, index);
      if (DustFluid) {
	WriteDiskPolar (dustvr, index);
	WriteDiskPolar (dustvt, index);
      }
    }
    //if (AdvecteLabel == YES) WriteDiskPolar (label, index);
    if (Write_Energy == YES) WriteDiskPolar (gasenerg, index);
    if (Write_Temperature == YES) WriteDiskPolar (Temperature, index);
    if (Write_DivV == YES) WriteDiskPolar (DivergenceVelocity, index);
    if (Write_ViscHeat == YES)  WriteDiskPolar (ViscHeat, index);
    if (Write_TherDiff == YES)  WriteDiskPolar (EntropyDiff, index);
    if (Write_RadDiff == YES)  WriteDiskPolar (RadiativeKCoeff, index);
    if (Write_TherCool == YES)  WriteDiskPolar (ThermCool, index);
    if (Write_Potential == YES)  WriteDiskPolar (Potential, index);
    if (Write_Test == YES)  WriteDiskPolar (Test, index);
    if (Write_gr == YES)  WriteDiskPolar (gr, index);
    if (Write_gtheta == YES) {
      WriteDiskPolar (gtheta, index);
      WriteDiskPolar (torquesg, index);
      WriteDiskPolar (torquesumdisc, index);
    }
    if ((Write_RadFBAcc == YES) && DustFeedback) WriteDiskPolar (RadFBAcc, index);
    if ((Write_AziFBAcc == YES) && DustFeedback) WriteDiskPolar (AziFBAcc, index);
    if ((Write_DustDensity == YES) && (NBPART != 0) ) WriteDiskPolar (dustpcdens, index);
    if ((Write_StokesNumber == YES) && (DustFluid) ) WriteDiskPolar (Stokes, index);
    MPI_Barrier (MPI_COMM_WORLD);
    if (Merge && (CPU_Number > 1)) merge (index);
    if (DustFluid && Merge && (CPU_Number > 1)) mergedust (index);
  }
}



void EmptyDustSystemFile(sys,timestep)
     DustSystem* sys;
     int timestep;
{
  FILE *output;
  char name[256];
  
  if (!CPU_Master) return;
  sprintf(name, "%sdustsystat%d.dat",OUTPUTDIR,timestep);
  output = fopen(name, "w");
  if (output == NULL) {
    fprintf (stderr, "Can't write dustsys file. Aborting \n");
    prs_exit (1);
  }
  fclose (output);
}


void WriteDustSystemFile (sys,timestep)
     DustSystem* sys;
     int timestep;
{
  FILE *output;
  char name[256];
  real r;
  extern boolean SelfGravity, SGZeroMode;
  char command[1024];
  int i, l, one_if_odd;
  masterprint("Updating dustsystat%d.dat...", timestep);
  fflush (stdout);
  if (CPU_Master)
    sprintf(name, "%sdustsystat%d.dat",OUTPUTDIR,timestep);
  else
    sprintf(name, "%sdustsystat%d.dat.%05d",OUTPUTDIR,timestep, CPU_Rank);
  output = fopen (name, "w");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'dustsystat%d.dat.%05d'. Aborting. \n",timestep, CPU_Rank);
    prs_exit (1);
  }
  for (i = 0; i<  NBPART; i++) {
    r = sys->r[i];
    if ( (r >= Rinf[Zero_or_active]) && (r < Rsup[Max_or_active-1]) ) {
      fprintf (output, "%#.6g\t%#.6g\t%#.6g\t%#.6g\t%#.6g\t%#.6g\n",sys->r[i],sys->th[i],sys->vr[i],sys->vth[i],sys->stokesnb[i],sys->dustsize[i]*unit_length);
    } 
    /*
    if ( (r < RMIN) || (r > RMAX) ) 
      printf ("Particle %d not written with r=%lg, th=%lg, vr=%lg, vth=%lg\n",i,sys->r[i],sys->th[i],sys->vr[i],sys->vth[i]);
    */
  }
  fclose (output);
  if (!CPU_Master) return;
  
  // Merging dustsyst files
  if ( SelfGravity && !SGZeroMode ) {
    one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
    for (i = 0; i < (CPU_Number+one_if_odd)/2; i++) {
      if ( i != 0 ) {
	sprintf (command, "cd %s; sync; cat dustsystat%d.dat.%05d >> dustsystat%d.dat", \
		 OUTPUTDIR, timestep, i, timestep);
	system (command);
      }
      l = (i + (CPU_Number + one_if_odd)/2)%(CPU_Number + one_if_odd);
      if ( i != CPU_Highest ) {
	sprintf (command, "cd %s; sync; cat dustsystat%d.dat.%05d >> dustsystat%d.dat", \
		 OUTPUTDIR, timestep, l, timestep);
	system (command);
      }
    }
    sprintf (command, "cd %s; rm -f dustsystat%d.dat.0*",	\
	     OUTPUTDIR, timestep);
    system (command);
  }
  else {
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; sync; cat dustsystat%d.dat.%05d >> dustsystat%d.dat", \
		 OUTPUTDIR, timestep, i, timestep);
      system (command);
    }
    sprintf (command, "cd %s; rm -f dustsystat%d.dat.0*",	\
	     OUTPUTDIR, timestep);
    system (command);
  }
  
  printf ("done\n");
  fflush(stdout);
}

void CreateJacobiFile(dsys,num)
     DustSystem* dsys;
     int num;
{
  FILE *output;
  char name[256];
  
  if (!CPU_Master) return;
  sprintf(name, "%sjacobi%d.dat",OUTPUTDIR,num);
  output = fopen(name, "w");
  if (output == NULL) {
    fprintf (stderr, "Can't write %d jacobi file. Aborting \n",num);
    prs_exit (1);
  }
  fclose (output);
}

void WriteJacobi(dsys, Plsys,  num)
     DustSystem *dsys;
     int num;
     PlanetarySystem *Plsys;
{
  FILE *output;
  char name[256];
  real xp, yp, rp, mp, tp, deltatheta;
  real rd, td, eps, d, PotPlan, OmegaPlan;
  fflush (stdout);
  sprintf(name, "%sjacobi%d.dat",OUTPUTDIR,num);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'jacobi%d.dat'. Aborting. \n",num);
    prs_exit (1);
  }
  xp = Plsys->x[0];
  yp = Plsys->y[0];
  tp = atan2(yp,xp);
  deltatheta = dsys->th[num] - tp; // theta_particle - theta_planet

  if (deltatheta < -M_PI)
    deltatheta += (2.0*M_PI);
  if (deltatheta > M_PI)
    deltatheta -= (2.0*M_PI);

  if (PhysicalTime == 0) {
    rd = dsys->r[num];
    td = dsys->th[num];
    rp = sqrt( xp*xp + yp*yp);
    mp = Plsys->mass[0];
    eps = compute_smoothing(rp);  // smoothing length at planet's orbital radius
    // Smoothed distance between planet and dust particle
    d = sqrt( rd*rd + rp*rp - 2.0*rd*rp*cos(td-tp) + eps*eps );
    PotPlan = mp/d;
    OmegaPlan = pow(rp,-1.5); 
    dsys->jacobi[num] = (0.5*(pow(dsys->vr[num],2.) + pow(dsys->vth[num],2.))) - 1./dsys->r[num] - PotPlan - OmegaPlan*dsys->r[num]*dsys->vth[num];
  }
  fprintf (output, "%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\t%#.14g\n",PhysicalTime,dsys->jacobi[num],dsys->r[num],deltatheta,dsys->th[num],dsys->vr[num],dsys->vth[num]);
  
  fclose (output);
  fflush(stdout);
}
