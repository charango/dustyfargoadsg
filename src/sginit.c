#include "mp.h"

void init_sg ()
{
  extern boolean LogGrid, SGZeroMode;
  if ( !LogGrid )
    erreur("A logarithmic grid is required to compute self-gravity with polar method. Start again!");
  SG_initcounter = 0;
  SGP_eps = SGTHICKNESSSMOOTHING*ASPECTRATIO;
  SGP_rstep = log(Radii[GLOBALNRAD]/Radii[0]) / (real)GLOBALNRAD;
  SGP_tstep = (PMAX-PMIN)/(real)NSEC;
  
  /* Fully self-gravitating case: 2D FFT are performed */
  if ( !SGZeroMode ) {
    SGP_Sr = (fftw_real *) fftw_malloc (sizeof(fftw_real) * total_local_size);
    SGP_St = (fftw_real *) fftw_malloc (sizeof(fftw_real) * total_local_size);
    SGP_Kr = (fftw_real *) fftw_malloc (sizeof(fftw_real) * total_local_size);
    SGP_Kt = (fftw_real *) fftw_malloc (sizeof(fftw_real) * total_local_size);
    SGP_Accr = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * total_local_size);
    SGP_Acct = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * total_local_size);
    if ( (SGP_Sr == NULL) || (SGP_St == NULL) || (SGP_Kr == NULL) || (SGP_Kt == NULL) ||  (SGP_Accr == NULL) || (SGP_Acct == NULL) ) {
      fprintf (stderr, "Not enough memory for allocation of SGP_tabs in sginit.c \n");
      prs_exit (1);
    }
    SG_Accr = (real *) malloc(sizeof(real) * hydro_totalsize);
    SG_Acct = (real *) malloc(sizeof(real) * hydro_totalsize);
    if ( (SG_Accr == NULL) || (SG_Acct == NULL) ) {
      fprintf (stderr, "Not enough memory for allocation of SG_Acc tabs in sginit.c \n");
      prs_exit (1);
    }
    /* Polargrids */
    gr = CreatePolarGrid(NRAD, NSEC, "sgaccr");
    gtheta = CreatePolarGrid(NRAD, NSEC, "sgacctheta");
    torquesg = CreatePolarGrid(NRAD, NSEC, "torquesg");
    /* Note that the creation of the forward fft plan is called in the
       d.d., in split.c */
    SGP_fftplan_backward = rfftw2d_mpi_create_plan(MPI_COMM_WORLD,
						   2*GLOBALNRAD, NSEC, 
						   FFTW_COMPLEX_TO_REAL, 
						   FFTW_MEASURE);
  }
  /* Axisymmetric component of the disk self-gravity only. The
     calculation is not parallelized, since only 1D FFT are
     required */
  if ( SGZeroMode ) {
    SGP_Sr = (fftw_real *) fftw_malloc (sizeof(fftw_real) * 2 * GLOBALNRAD);
    SGP_Kr = (fftw_real *) fftw_malloc (sizeof(fftw_real) * 2 * GLOBALNRAD);
    SGP_Accr = (fftw_real *) fftw_malloc (sizeof(fftw_real) * 2 * GLOBALNRAD);
    SG_Accr = (real *) malloc(sizeof(real) * GLOBALNRAD);
    if ( (SGP_Sr == NULL) || (SGP_Kr == NULL) || (SGP_Accr == NULL) ||  (SG_Accr == NULL) ) {
      fprintf (stderr, "Not enough memory for allocation of SGZeroMode tabs in sginit.c \n");
      prs_exit (1);
    }
    GLOBAL_Axidens = (real *) malloc(sizeof(real) * GLOBALNRAD);
    if ( GLOBAL_Axidens==NULL ) {
      fprintf (stderr, "Not enough memory for allocation of GLOBAL_Axidens in sginit.c \n");
      prs_exit (1);
    }
    SGP_fft1Dplan_forward  = rfftw_create_plan(2*GLOBALNRAD,
					       FFTW_REAL_TO_COMPLEX,
					       FFTW_MEASURE|FFTW_IN_PLACE);
    SGP_fft1Dplan_backward = rfftw_create_plan(2*GLOBALNRAD,
					       FFTW_COMPLEX_TO_REAL,
					       FFTW_MEASURE|FFTW_IN_PLACE);
  }
  GLOBAL_AxiSGAccr = (real *) malloc(sizeof(real) * GLOBALNRAD);
  if ( GLOBAL_AxiSGAccr == NULL ) {
    fprintf (stderr, "Not enough memory for allocation of GLOBAL_AxiSGAccr in sginit.c \n");
    prs_exit (1);
  }
}
