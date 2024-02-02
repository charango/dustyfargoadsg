#include "mp.h"

void compute_SGZeroMode (Rho)
     PolarGrid *Rho;
{
  int i, j, mr, mi;
  real u, normaccr;
  real num, den, integrale;
  real *dens;
  dens = Rho->Field;
  mpi_make1Dprofile (dens, GLOBAL_Axidens);
  
  /* Reduced density and its 1D Fourier transform. */
  for ( i = 0; i < GLOBALNRAD; i++ )
    SGP_Sr[i] = GLOBAL_Axidens[i] * sqrt (GlobalRmed[i] / GlobalRmed[0]);
  for ( i = GLOBALNRAD; i < 2*GLOBALNRAD; i++ )
    SGP_Sr[i] = 0.;
  rfftw_one (SGP_fft1Dplan_forward, SGP_Sr, NULL);
  SGP_buffft_Sr = (real *) SGP_Sr;
  
  /* Reduced kernel and its 1D Fourier transform. */
  if ( SG_initcounter == 0 ) {
    for ( i = 0; i < 2*GLOBALNRAD; i++ ) {
      if ( i < GLOBALNRAD )
	u = log (Radii[i] / Radii[0]);
      else
	u = -log(Radii[2*GLOBALNRAD-i]/Radii[0]);
      SGP_Kr[i] = 0.;
      for ( j = 0; j < NSEC; j++ ) {
	SGP_Kr[i] += ( 1.0 + SGP_eps*SGP_eps - CosAzimuth[j]*exp(-u) ) * \
	  pow(SGP_eps*SGP_eps*exp(u) + 2.0*( cosh(u) - CosAzimuth[j] ),-1.5);
      }
      SGP_Kr[i] *= SGP_tstep;
    }
    rfftw_one (SGP_fft1Dplan_forward, SGP_Kr, NULL); 
    SGP_buffft_Kr = (real *) SGP_Kr;
  }
  
  /* Convolution product of TF(Sr) and TF(Kr) in Fourier space. Be 
     careful that in 1D, contrary to 2D case, complex elements 
     are ordered in 'halfcomplex' way, with all real parts then ima-
     ginary parts. Here, the total array dimension is 2*GLOBALNRAD 
     thus even, so the order is Re(X0), Re(X1), ..., Re(Xn), Im(Xn-1),
     ..., Im(X1). See fftw online manual for more details. */
  SGP_buffft_Accr = (real *) SGP_Accr;
  SGP_buffft_Accr[0] = -G*( SGP_buffft_Kr[0]*SGP_buffft_Sr[0] );
  for ( i = 1; i < GLOBALNRAD; i++ ) { // real parts
    mr = i;
    mi = 2*GLOBALNRAD - i;
    SGP_buffft_Accr[i] = -G*( SGP_buffft_Kr[mr]*SGP_buffft_Sr[mr] -	\
			      SGP_buffft_Kr[mi]*SGP_buffft_Sr[mi] );
  }
  mr = GLOBALNRAD;
  SGP_buffft_Accr[mr] = -G*( SGP_buffft_Kr[mr]*SGP_buffft_Sr[mr] );
  for ( i = GLOBALNRAD + 1; i < 2*GLOBALNRAD; i++ ) { // imaginary parts
    mr = 2*GLOBALNRAD - i;
    mi = i;
    SGP_buffft_Accr[i] = -G*( SGP_buffft_Kr[mr]*SGP_buffft_Sr[mi] +	\
			      SGP_buffft_Kr[mi]*SGP_buffft_Sr[mr] );
  }
  
  /* TF-1 of TF(Accr). */
  rfftw_one (SGP_fft1Dplan_backward, SGP_Accr, NULL);
  SGP_buffft_Accr = (real *) SGP_Accr;
  
  /* Normalization. */
  for ( i = 0; i < GLOBALNRAD; i++ ) {
    SG_Accr[i] = SGP_buffft_Accr[i];
    normaccr = SGP_rstep / ( 2.0*(real)GLOBALNRAD );
    normaccr /= sqrt(GlobalRmed[i] / GlobalRmed[0]);
    SG_Accr[i] *= normaccr;
  }
  
  /* Compensation from self-forces. */
  integrale = 0.0;
  for ( j = 0; j < NSEC; j++) {
    num = 1.0 + SGP_eps*SGP_eps - CosAzimuth[j];
    den = pow(2.0*(1.0-CosAzimuth[j]) + SGP_eps*SGP_eps,1.5);
    integrale += num/den;
  }
  integrale *= SGP_tstep;
  for ( i = 0; i < GLOBALNRAD; i++ )
    SG_Accr[i] += G*GLOBAL_Axidens[i]*SGP_rstep*integrale;
}
