#include "mp.h"

void compute_selfgravity (Rho, RadialVelocity, AzimutalVelocity, DustRadialVelocity, DustAzimutalVelocity, DeltaT, SGUpdate)
     PolarGrid *Rho;
     PolarGrid *RadialVelocity;
     PolarGrid *AzimutalVelocity;
     PolarGrid *DustRadialVelocity, *DustAzimutalVelocity;
     real DeltaT;
     boolean SGUpdate;
{
  extern boolean SGZeroMode;
  if ( SG_initcounter == 0 ) {
    init_sg ();
    if ( !SGZeroMode )
      compute_fftkernel ();
  }
  /* Only the axisymmetric component of the disk self-gravity is taken
     into account */
  if ( SGZeroMode )
    compute_SGZeroMode (Rho);
  /* Fully self-gravitating case */
  if ( !SGZeroMode ) {
    /* Update calculation of FFT (reduced density) */
    compute_fftdensity (Rho);
    /* We now compute radial and azimutal components of sg
       acceleration as a convolution product of reduced density and
       kernel arrays */
    compute_sgacc (Rho);
  }
  if ( SGUpdate ) {
    /* Updates values of vrad and vtheta with self-gravity */
    update_sgvelocity (RadialVelocity, AzimutalVelocity, DustRadialVelocity, DustAzimutalVelocity, DeltaT);
  }
  SG_initcounter++;
}
