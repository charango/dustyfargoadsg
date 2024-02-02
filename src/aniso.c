#include "mp.h"
/* In this function we calculate the anisotropy coefficient alpha that
   reproduces the non-axisymmetric component of the disk self-gravity,
   by means of an anisotropic pressure tensor (see Baruteau & Masset
   08b). Alpha is defined here as alpha = 1-beta/Q, where beta depends
   only on the smoothing length to disk thickness ratio
   (eta=eps/H): */
/*           eta           beta       */
/*           0.1          0.32(4)     */
/*           0.3          0.61(4)     */
/*           0.6          0.94(1)     */

void compute_anisotropic_pressurecoeff (sys)
     PlanetarySystem *sys;
{
  real xpl, ypl, rpl, mpl;
  real Q, beta, alpha;
  mpl = sys->mass[0];
  xpl = sys->x[0];
  ypl = sys->y[0];
  rpl = sqrt(xpl*xpl + ypl*ypl);
  Q = pow(rpl,-2.0+SIGMASLOPE)*ASPECTRATIO/M_PI/SIGMA0;
  beta = 0.0;
  if ( fabs(SGTHICKNESSSMOOTHING-0.1) < 1e-2 ) beta=0.324;
  if ( fabs(SGTHICKNESSSMOOTHING-0.3) < 1e-2 ) beta=0.614;
  if ( fabs(SGTHICKNESSSMOOTHING-0.6) < 1e-2 ) beta=0.941;
  alpha = 1.0 - beta/Q;
  if (SG_initcounter == 1) {
    printf ("Q = %lg, beta = %lg et alpha = %lg\n", Q, beta, alpha);
  }
  SG_aniso_coeff = alpha;
}
