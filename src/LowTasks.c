/** \file LowTasks.c

Contains many low level short functions.  The name of these functions
should be self-explanatory in most cases.  The prefix 'prs_' stands
for 'personal'. The prefix 'master' means that only the process 0
executes the function [note that the architecture is not of the kind
master/slaves, all processes perform similar tasks, but a minor number
of tasks (like output of information on the standard output) do not
need to be performed by several processes.] The function fopenp() is
an upper layer of fopen(), which should be used only in the case of
writing or appending a file (and not reading a file). It tries to
create the output directory if it does not exist, and it issues an
error message if it fails, so that the calling function does not need
to worry about these details.
*/

#include "mp.h"
#include <stdarg.h>

real GetGlobalIFrac (r)
     real r;
{
  int i=0;
  real ifrac;
  if (r < GlobalRmed[0]) return 0.0;
  if (r > GlobalRmed[GLOBALNRAD-1]) return (real)GLOBALNRAD-1.0;
  while (GlobalRmed[i] <= r) i++;
  ifrac = (real)i+(r-GlobalRmed[i-1])/(GlobalRmed[i]-GlobalRmed[i-1])-1.0;
  return ifrac;
}

void prs_exit (numb)
     int numb;
{
  MPI_Finalize ();
  exit (numb);
}

void masterprint (const char *template, ...)
{
  va_list ap;
  if (!CPU_Master) return;
  va_start (ap, template);
  vfprintf (stdout, template, ap);
  va_end (ap);
}

void mastererr (const char *template, ...)
{
  va_list ap;
  if (!CPU_Master) return;
  va_start (ap, template);
  vfprintf (stderr, template, ap);
  va_end (ap);
}

void erreur(string)
     char *string;
{
  fprintf(stderr, "%s\n", string);
  prs_exit(1);
}

void *prs_malloc (number_of_bytes)
     size_t number_of_bytes;
{
  void *ptr;
  long i;
  if (number_of_bytes <= 0)
    return NULL;
  ptr = malloc (number_of_bytes);
  if (ptr == NULL) erreur ("Not enough memory.");
  for (i = 0; i < number_of_bytes; i++)
    *((char *)ptr+i) = 0;
  return ptr;
}

void message (msg) 
     char *msg;
{
  fprintf (stdout, "%s", msg);
}

PolarGrid    *
CreatePolarGrid(Nr, Ns, name)
int             Nr, Ns;
char           *name;
{
  PolarGrid      *array;
  real           *field;
  char           *string;
  int             i, j, l;
  
  array = (PolarGrid *) malloc(sizeof(PolarGrid));
  if (array == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  field = (real *) malloc(sizeof(real) * (Nr + 3) * (Ns + 1) + 5);
  if (field == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  string = (char *) malloc(sizeof(char) * 256);
  if (string == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  sprintf(string, "%s", name);
  array->Field = field;
  array->Name = string;
  array->Nrad = Nr;
  array->Nsec = Ns;
  for (i = 0; i <= Nr; i++) {
    for (j = 0; j < Ns; j++) {
      l = j + i*Ns;
      field[l] = 0.;
    }
  }
  return array;
}


void MultiplyPolarGridbyConstant (arraysrc, constant)
     PolarGrid *arraysrc;
     real constant;
{
  int i, nr, ns;
  real *fieldsrc;
  nr = arraysrc->Nrad;
  ns = arraysrc->Nsec;
  fieldsrc  =  arraysrc->Field;
#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) {
    fieldsrc[i] *= constant;
  }
}

void DumpSources (argc, argv, parfile, cfgfile)
     int argc;
     char *argv[], *parfile, *cfgfile;
{
  char CommandLine[1024];
  char filecom[1024];
  int i;
  FILE *COM;
  if (!CPU_Master) return;
  sprintf (CommandLine, "cp source.tar.bz2 %ssrc.tar.bz2", OUTPUTDIR);
  system (CommandLine);
  sprintf (CommandLine, "cp %s %s %s", parfile, cfgfile, OUTPUTDIR);
  system (CommandLine);
  sprintf (filecom, "%srun.commandline", OUTPUTDIR);
  COM = fopen (filecom, "w");
  if (COM == NULL) {
    mastererr ("Could not open %s\nAborted.\n", filecom);
    prs_exit(1);
  }
  for (i = 0; i < argc; i++) {
    fprintf (COM, "%s ",argv[i]);
  }
  fclose (COM);
}


void MakeDir (string)
     char *string;
{
  int foo=0;
  char command[MAX1D];
  DIR *dir;
  /* Each processor tries to create the directory, sequentially */
  /* Silent if directory exists */
  if (CPU_Rank) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 53, MPI_COMM_WORLD, &fargostat);
  dir = opendir (string);
  if (dir) {
    closedir (dir);
  } else {
    fprintf (stdout, "Process %d creates the directory %s\n", CPU_Rank, string);
    sprintf (command, "mkdir -p %s", string);
    system (command);
  }
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 53, MPI_COMM_WORLD);
}


FILE *fopenp (string, mode)
     char *string, *mode;
{
  FILE *f;
  f = fopen (string, mode);
  if (f == NULL) {
    /* This should be redundant with the call to MakeDir () at the
       beginning, from main.c; this is not a problem however */
    printf ("Process %d could not open %s\n", CPU_Rank, string);
    printf ("Trying to create %s\n", OUTPUTDIR);
    MakeDir (OUTPUTDIR);
    f = fopen (string, "w");	/* "w" instead of mode: at this stage we know the file does not exist */
    if (f == NULL) {
      fprintf (stdout, "I still cannot open %s.\n", string);
      fprintf (stdout, "You should check that the permissions are correctly set.\n");
      fprintf (stdout, "Run aborted\n");
      prs_exit (1);
    }
  }
  return f;
}

void ComputeCodeUnits ()
{
  char filename[256];
  FILE 	*dim;
  /* by default, unit_mass is one Sun mass in kg */
  unit_mass   = 1.9891e30*FACTORUNITMASS;     
  /* by default, unit_length is one AU in meters */
  unit_length = 1.49598e11*FACTORUNITLENGTH;  
  /* by default, mmw is 2.35 (Solar System composition) */
  mmw         = 2.35*FACTORMMW;              
  unit_time = sqrt( pow(unit_length,3.) / 6.673e-11 / unit_mass );
  unit_temperature = mmw * 8.0841643e-15 * unit_mass / unit_length;
  // This is Stefan-Boltzmann constant
  sigma_SB = 5.6704e-8 * pow(unit_mass,-1) * pow(unit_time,3.) * pow(unit_temperature,4.);
  //masterprint ("sigma_SB = %lg",sigma_SB);
  masterprint("Length unit:      %lg\n",unit_length);
  masterprint("Mass unit:        %lg\n",unit_mass);
  masterprint("Time unit:        %lg\n",unit_time);
  masterprint("Temperature unit: %lg\n",unit_temperature);
  if (!CPU_Master) return;
  sprintf (filename, "%sunits.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    prs_exit (1);
  }
  fprintf (dim,"%.18g\t%.18g\t%.18g\t%.18g\n",unit_mass,unit_length,unit_time,unit_temperature);
  fclose (dim);
}
