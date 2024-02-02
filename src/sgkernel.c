#include "mp.h"

void compute_fftkernel ()
{
  int i, j, l;
  int stride;
  real u;
  stride = 2*(NSEC/2+1);
  
  for ( i = 0; i < local_Nx; i++ ) {
    if ( i+local_i_start < GLOBALNRAD )
      u = log(Radii[i+local_i_start]/Radii[0]);
    else
      u = -log(Radii[2*GLOBALNRAD-(i+local_i_start)]/Radii[0]);
    for ( j = 0; j < NSEC; j++ ) {
      l = i*stride + j;
      SGP_Kr[l] = 1.0 + SGP_eps*SGP_eps - CosAzimuth[j]*exp(-u);
      SGP_Kr[l] *= pow(SGP_eps*SGP_eps*exp(u) + 2.0*( cosh(u) - CosAzimuth[j] ),-1.5);
      SGP_Kt[l] = SinAzimuth[j];
      SGP_Kt[l] *= pow(SGP_eps*SGP_eps*exp(u) + 2.0*( cosh(u) - CosAzimuth[j] ),-1.5);
    }
  }
  rfftwnd_mpi (SGP_fftplan_forward, 1, SGP_Kr, NULL, FFTW_TRANSPOSED_ORDER);
  SGP_buffft_Kr = (real *) SGP_Kr;
  rfftwnd_mpi (SGP_fftplan_forward, 1, SGP_Kt, NULL, FFTW_TRANSPOSED_ORDER);
  SGP_buffft_Kt = (real *) SGP_Kt;
}
