#include "mp.h"

void update_sgvelocity (VRad, VTheta, DVRad, DVTheta, Dt)
     PolarGrid *VRad, *VTheta;
     PolarGrid *DVRad, *DVTheta;
     real Dt;
{
  extern boolean SGZeroMode, DustFluid, DustFeelSG, DustFeelSGZeroMode;
  int i, j, l, nr;
  int jm1, lm1;
  real *vrad, *vtheta;
  real *dvrad, *dvtheta;
  real *radsgacc, *azisgacc;
  vrad = VRad->Field;
  vtheta = VTheta->Field;
  nr = VTheta->Nrad;
  radsgacc = RadSGAcc->Field;
  azisgacc = AziSGAcc->Field;
  if (DustFluid) {
    dvrad = DVRad->Field;
    dvtheta = DVTheta->Field;
  }
  
  /* Particular case where particles only feel the azisymmetric part
     of the disc's self-gravitating potential */
  if (DustFluid && DustFeelSGZeroMode)
    mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);

  /* Here we update velocity fields with self-gravity */
  for ( i = 0 ; i < nr; i++ ) {
    for (j = 0; j < NSEC; j++) {
      l = i*NSEC + j;
      /* We compute VRAD - half-centered in azimuth - from
	 centered-in-cell radial sg acceleration, as a linear
	 interpolation */
      if ( i > 0 ) {
	if ( !SGZeroMode ) {
	  radsgacc[l] = ( (Rinf[i] - Rmed[i-1]) * SG_Accr[l] +		\
			  (Rmed[i] - Rinf[i]) * SG_Accr[l-NSEC] ) * InvDiffRmed[i];
	  vrad[l] += Dt*radsgacc[l];
	  if (DustFluid && DustFeelSG && (!DustFeelSGZeroMode)) 
	    dvrad[l] += Dt*radsgacc[l];
	  if (DustFluid && DustFeelSGZeroMode) {
	    radsgacc[l] = ( (Rinf[i] - Rmed[i-1]) * GLOBAL_AxiSGAccr[i+IMIN] + \
			    (Rmed[i] - Rinf[i]) * GLOBAL_AxiSGAccr[i+IMIN-1] ) * InvDiffRmed[i];
	    dvrad[l] += Dt*radsgacc[l];
	  }
	}
	if ( SGZeroMode ) {
	  radsgacc[l] = ( (Rinf[i] - Rmed[i-1]) * SG_Accr[i+IMIN] +	\
			  (Rmed[i] - Rinf[i]) * SG_Accr[i+IMIN-1] ) * InvDiffRmed[i];
	  vrad[l] += Dt*radsgacc[l];
	  if (DustFluid && DustFeelSG) 
	    dvrad[l] += Dt*radsgacc[l];
	}
      }
      /* We compute VTHETA - half-centered in radius - from
	 centered-in-cell azimutal sg acceleration, as a linear
	 interpolation */
      if ( !SGZeroMode ) {
	if (j == 0) 
	  jm1 = NSEC-1;
	else
	  jm1 = j-1;
	lm1 = i*NSEC + jm1;
	azisgacc[l] = 0.5*(SG_Acct[l] + SG_Acct[lm1]);
	vtheta[l] += Dt*azisgacc[l];
	if (DustFluid && DustFeelSG && (!DustFeelSGZeroMode)) 
	  dvtheta[l] += Dt*azisgacc[l];
      }
    }
  }
}
