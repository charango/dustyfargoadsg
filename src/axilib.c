#include "mp.h"

/* This function calculates a global, axisymmetric field ("axifield")
   from current field "gridfield" */
void mpi_make1Dprofile (gridfield, axifield)
     real *gridfield;
     real *axifield;
{
  MPI_Request req1;
  int i, j, l;
  real *localaxifield;
  localaxifield = (real*) malloc(sizeof(real) * NRAD);
  if ( localaxifield == NULL ) 
    erreur ("Not enough memory in axilib.c ; suspicious...");
  
  /* We first calculate axisymmetric local field */
  for ( i = 0; i < NRAD; i++ )
    localaxifield[i] = 0.;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      localaxifield[i] += gridfield[l];
    }
    localaxifield[i] /= (real)NSEC;
  }
  /* Then we share it with other cpus to yield a global, axisymmetric
     field */
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD; i++ )
      axifield[i] = localaxifield[i];
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for ( i = 0; i < GLOBALNRAD; i++ ) {
	if ( i < Max_or_active )
	  axifield[i] = localaxifield[i];
	else
	  axifield[i] = 0.;
      }
      MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (i = Zero_or_active; i < Max_or_active; i++)
	axifield[i+IMIN] = localaxifield[i];
      if ( CPU_Rank != CPU_Highest ) {
	MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait (&req1, &fargostat);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  }
  free (localaxifield);
}


/* This function calculates a global field ("globalfield") from
   current field "gridfield" */
void mpi_makeglobalfield (gridfield, globalfield)
     real *gridfield;
     real *globalfield;
{
  MPI_Request req1;
  int i, j, l;
  real *localglobfield;
  localglobfield = (real*) malloc(sizeof(real) * NRAD * NSEC);
  if ( localglobfield == NULL ) 
    erreur ("Not enough memory for localglobfield in axilib.c ; suspicious...");
  
  /* We first calculate local "global field" */
  for ( i = 0; i < NRAD*NSEC; i++ )
    localglobfield[i] = 0.;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      localglobfield[l] = gridfield[l];
    }
  }
  /* Then we share it with other cpus to yield a global field on the polar grid */
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD*NSEC; i++ )
      globalfield[i] = localglobfield[i];
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for (i = 0; i < GLOBALNRAD; i++) {
	for (j = 0; j < NSEC; j++) {
	  l = i*NSEC + j;
	  if ( i < Max_or_active )
	    globalfield[l] = localglobfield[l];
	  else
	    globalfield[l] = 0.;
	}
      }
      MPI_Isend (globalfield, GLOBALNRAD*NSEC, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (globalfield, GLOBALNRAD*NSEC, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (i = Zero_or_active; i < Max_or_active; i++) {
	for (j = 0; j < NSEC; j++) {
	  l = i*NSEC + j;
	  globalfield[(i+IMIN)*NSEC+j] = localglobfield[l];
	}
      }
      if ( CPU_Rank != CPU_Highest ) {
	MPI_Isend (globalfield, GLOBALNRAD*NSEC, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait (&req1, &fargostat);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (globalfield, GLOBALNRAD*NSEC, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  }
  free (localglobfield);
}
