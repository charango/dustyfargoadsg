#include "mp.h"

/* Static variables for turbulent potential expression */
static int m[MAX1D];
static real xi[MAX1D], rc[MAX1D], sigma[MAX1D], phic[MAX1D], omegac[MAX1D], ttilde[MAX1D], t0[MAX1D], tf[MAX1D], deltat[MAX1D];

/* The routine below adds wave-like stochastic modes to the
gravitational potential felt by the disk, with statistical properties
aimed at reproducing those of actual disk turbulence driven by the
MRI. It uses the method originally developed by Laughlin et al. (2004,
LSA04) and modified by Baruteau & Lin (2010, BL10). */

void ApplyLSAOnPotential ()
{
  int i, j, l, k, nr, ns;
  real ampl, x1, x2, w, gamma, r, angle, mygamma, rmin, rmax;
  real *Pot, *TurbPot;
  FILE *fich;
  extern boolean HighMCutoff;
  char filename[200];
  Pot = Potential->Field;
  TurbPot = TurbPotential->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  for (i = 0; i < (nr+1)*ns; i++) TurbPot[i] = 0.0;
  /* The expression of the turbulent potential can be found, e.g., in
     Eqs 1 and 2 in BL10. Note that only CPU_Highest randomly sorts
     rc, m, phic and xi */
  if (CPU_Rank == CPU_Highest) {
    for (k = 0; k < NBTURBMODES; k++) {
      if (PhysicalTime >= tf[k]) {
	/* We sort the logarithm of the wavenumber m uniformly between
	   log(1) and log(Nsec/8) */

	// CB: improve by specifying m_min and m_max. m_min should be
	// set by the grid's azimuthal extent PMAX-PMIN, and m_max by
	// the minimum length scale that should be resolved, like for
	// instance H.
      
	m[k] = (long)exp(log((real)ns/8.0)*drand48());
	/* We sort rc uniformly throughout the grid except in the
	 wave-killing zones; TURBRMIN AND TURBRMAX are the radial 
	 boundaries where stochastic forcing is applied in the disc. */
	rc[k] = TURBRMIN + (TURBRMAX-TURBRMIN)*drand48();
	/* PMIN AND PMAX the min and max azimuths in the computational
	   grid (default: 0 and 2pi) */
	phic[k] = (PMAX-PMIN)*drand48();
	/* Parameter xi is sorted with a gaussian distribution of mean
	   0 and standard deviation unity */
	do {
	  x1 = 2.0*drand48()-1.0;
	  x2 = 2.0*drand48()-1.0;
	  w = x1*x1 + x2*x2;
	} while ( w >= 1.0);
	w = sqrt( (-2.0 * log(w)) / w );
	xi[k] = x1*w;
      }
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  /* CPU_Highest broadcasts information of turbulent parameters to all
     other cpus */
  MPI_Bcast (m, NBTURBMODES, MPI_INT, CPU_Highest, MPI_COMM_WORLD);
  MPI_Bcast (rc, NBTURBMODES, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  MPI_Bcast (phic, NBTURBMODES, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  MPI_Bcast (xi, NBTURBMODES, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  /* All cpus now infer the other modes' properties (simple algebra) */
  /* case where a cavity is implemeted (see below, default: rmin=0,
     rmax=H) */
  rmin = CAVITYRADIUS-CAVITYWIDTH*ASPECTRATIO;   
  rmax = CAVITYRADIUS+CAVITYWIDTH*ASPECTRATIO;
  for (k = 0; k < NBTURBMODES; k++) {
    if (PhysicalTime >= tf[k]) {
      omegac[k] = pow(rc[k],-1.5);
      sigma[k] = M_PI*rc[k]/4.0/(real)m[k];
      /* This expression of deltat uses the isothermal sound speed */
      deltat[k] = LSAMODESPEEDUP*2.0*M_PI*rc[k]/(real)m[k]/ASPECTRATIO/pow(rc[k],-0.5+FLARINGINDEX);
      /* Modes lifetime multiplied by a factor LSAMODESPEEDUP compared with LSA 04 */
      t0[k] = PhysicalTime;         // time at which a mode becomes active
      tf[k] = t0[k] + deltat[k];    // time at which a mode becomes inactive
      /* 
      // should be written at same frequency as gas*.dat files also,
      // the restart case should be worked out
      if (CPU_Rank == CPU_Highest) {
	sprintf (filename, "%sturb.dat", OUTPUTDIR);
	fich = fopen(filename, "a");
	fprintf(fich,"%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", k, m[k], rc[k], deltat[k], t0[k], tf[k], phic[k], xi[k], PhysicalTime);
	fclose(fich);
      }
      */
    }
    ttilde[k] = PhysicalTime-t0[k];
    /* We then update the grav. potential felt by the disc */
    gamma = GAMMATURB;
    /* when a cavity is accounted for, this is the minimum value of gamma */
    if ( (m[k] <= 6) || (!HighMCutoff) ) {
      /* If HighMCutoff is set to yes, modes with azimuthal wavenumber m>6 are discarded */
      for (i = 0; i < nr; i++) {
	r = Rmed[i];
	/* Heaviside function for gamma coefficient */
	if (r < rmin) mygamma = gamma*CAVITYRATIO;  // !default: cavityratio=1
	if ((r >= rmin) && (r <= rmax)) {
	  mygamma = gamma*exp((rmax-r)/(rmax-rmin)*log(CAVITYRATIO));
	}
	if (r > rmax) mygamma = gamma;
	ampl = mygamma/r*xi[k]*					\
	  exp( -(r-rc[k])*(r-rc[k])/sigma[k]/sigma[k] )*	\
	  sin( M_PI*ttilde[k]/deltat[k] );
	for (j = 0; j < ns; j++) {
	  angle = Azimuth[j];
	  l = j+i*ns;
	  TurbPot[l] += ampl*cos( m[k]*angle-phic[k]-(omegac[k]-OmegaFrame)*ttilde[k] );
	  Pot[l] += ampl*cos( m[k]*angle-phic[k]-(omegac[k]-OmegaFrame)*ttilde[k] );
	}
      }
    }
  }
}
