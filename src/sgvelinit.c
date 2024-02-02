#include "mp.h"

void init_azimutalvelocity_withSG (Vtheta)
     PolarGrid *Vtheta;
{
  extern boolean SGZeroMode;
  int i, j, l, ns, nr;
  real r, invr;
  real omega, omegakep;
  real *vt;
  vt = Vtheta->Field;
  nr  = Vtheta->Nrad;
  ns  = Vtheta->Nsec;
  if ( !SGZeroMode )
    mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
  else
    GLOBAL_AxiSGAccr = SG_Accr;
  for (i = 0; i < nr; i++) {
    r = Rmed[i];
    invr = 1./Rmed[i];
    omegakep = sqrt(G*1.0*invr*invr*invr);
    omega = sqrt( omegakep*omegakep*( 1.0 -				\
				      (1.+SIGMASLOPE-2.0*FLARINGINDEX)*pow(ASPECTRATIO,2.0)*pow(r,2.0*FLARINGINDEX) ) - \
		  invr*GLOBAL_AxiSGAccr[i+IMIN] );
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      vt[l] = r*omega;
    }
  }
}
