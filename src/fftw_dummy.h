/***********************************************/
/*                                             */
/*                                             */
/*  Fake functions library for non fftw built  */
/*                                             */
/*                                             */
/***********************************************/
#ifndef __APPLE__
#include <malloc.h>
#endif

#define FFTW_TRANSPOSED_ORDER 1
#define FFTW_COMPLEX_TO_REAL  1
#define FFTW_MEASURE  1
#define FFTW_REAL_TO_COMPLEX 1
#define FFTW_IN_PLACE 1

typedef int rfftwnd_mpi_plan;
typedef int rfftw_plan;
typedef double fftw_real;
typedef double fftw_complex;

void * fftw_malloc();
double rfftw2d_mpi_create_plan();
void rfftwnd_mpi();
void rfftwnd_mpi_destroy_plan();
double rfftw_create_plan();
void rfftw_one();
void rfftwnd_mpi_local_sizes();
