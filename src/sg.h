/* ------------------------------------- */
/*       Global variables used for       */ 
/*   computation of disk self-gravity    */
/* ------------------------------------- */

int SG_initcounter;      /* Initialization counter */

/* ================================= */
/* Self-gravity on a polar grid */
/* ================================= */

real SGP_eps;             /* Smoothing parameter of potential */
real SGP_rstep;
real SGP_tstep;
real SG_aniso_coeff;   /* to mimic non-axisymmetric effects of self-gravity */
/* with an anisotropic pressure */

/* The plans for 2D and 1D fft computation */
rfftwnd_mpi_plan SGP_fftplan_forward;
rfftwnd_mpi_plan SGP_fftplan_backward;
rfftw_plan SGP_fft1Dplan_forward;
rfftw_plan SGP_fft1Dplan_backward;

/* Arrays and fft-arrays, since (mpi)fft computes in-place transforms... */
fftw_real *SGP_Sr, *SGP_St;
fftw_real *SGP_Kr, *SGP_Kt;
fftw_complex *SGP_Accr, *SGP_Acct;
real *SG_Accr, *SG_Acct;

/* Buffer arrays needed for the fft computation */
real *SGP_buffft_Sr, *SGP_buffft_St;
real *SGP_buffft_Kr, *SGP_buffft_Kt;
real *SGP_buffft_Accr, *SGP_buffft_Acct;

/* Global axisymetric arrays used in initialization or for zero mode of SG */
real *GLOBAL_Axidens;
real *GLOBAL_AxiSGAccr;
