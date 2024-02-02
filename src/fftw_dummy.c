/***********************************************/
/*                                             */
/*                                             */
/*  Fake functions library for non fftw built  */
/*                                             */
/*                                             */
/***********************************************/

#include "fftw_dummy.h"
#include <stdio.h>
#include <stdlib.h>

void * fftw_malloc(a)
     int a;
{
  double * b;
  b = (double*) malloc(a);
  return b ;
}

double rfftw2d_mpi_create_plan()
{
  return 1.;
}

void rfftwnd_mpi()
{
}

void rfftwnd_mpi_destroy_plan()
{
}

double rfftw_create_plan()
{
  return 1.;
}

void rfftw_one()
{
}

void rfftwnd_mpi_local_sizes()
{
}
