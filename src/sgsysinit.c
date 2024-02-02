#include "mp.h"

void init_planetarysys_withSG (sys)
     PlanetarySystem *sys;
{
  extern boolean SGZeroMode;
  int k, ipl;
  real x, y, r, dist, ri, rip1, dr, sgacc;
  
  /* Here we calculate global, axisymmetric density field, known by
     all cpus */
  if ( !SGZeroMode )
    mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
  else
    GLOBAL_AxiSGAccr = SG_Accr;
  
  /* Planetary system initialization with self-gravity: planets are
     put in a fixed circular orbit, we need to know radial sg
     acceleration felt by planets. */
  for ( k = 0; k < sys->nb; k++ ) {
    x = sys->x[k];
    y = sys->y[k];
    r = sqrt(x*x + y*y);
    /* dist denotes the planet's semi-major axis */
    dist = (real)(r / (1. + ECCENTRICITY));
    
    if (dist < GlobalRmed[GLOBALNRAD-1]) {
      ipl = 0;
      while ( (GlobalRmed[ipl] <= dist) && (ipl < GLOBALNRAD-2) ) ipl++;
      ri = GlobalRmed[ipl];
      rip1 = GlobalRmed[ipl+1];
      dr = rip1 - ri;
      sgacc = (dist - ri)*GLOBAL_AxiSGAccr[ipl+1] + (rip1 - dist)*GLOBAL_AxiSGAccr[ipl];
      sgacc /= dr;
    }
    else {
      sgacc = GLOBAL_AxiSGAccr[GLOBALNRAD-1];
    }
    
    /* sgacc is the radial sg acc. at the planet's semi-major axis */
    /* We disentangle between binary and non-binary cases */
    if ( (sys->Binary[k] == YES) && (k == 0) )
      sys->vy[k] *= (real)sqrt (1. - dist*dist*sgacc/G/(1.0+sys->mass[k]+sys->mass[k+1]) );
    if ( (sys->Binary[k] == YES) && (k == 1) )
      sys->vy[k] *= (real)sqrt (1. - dist*dist*sgacc/G/(1.0+sys->mass[k]+sys->mass[k-1]) );
    if (sys->Binary[k] == NO)
      sys->vy[k] *= (real)sqrt (1. - dist*dist*sgacc/G/(1.0+sys->mass[k]) );
  }
}
