/** \file Pframeforce.c

Functions that evaluate the %force between the planets and the disk.
The FillForcesArrays() function is ill-named: it rather fill an array
of the potential, that is later use to derive the force acting on the
disk at every zone.  The name of this file is due to the fact that we
work in the frame centered on the primary (which is therefore not
inertial). Old versions of fargo also featured the file Gframeforce.c,
in which we worked in the frame centered on the center of gravity of
the system.  The present file also contains the functions necessary to
update the planets' position and velocities (taking into account, or
not, the indirect term, ie the term arising from the fact that the
frame is not inertial), as well as the function that initializes the
hydrodynamics fields with analytic prescription.

*/

#include "mp.h"

extern boolean AllowAccretion, Indirect_Term, Discard_GasIndirect_term, DustFluid, DustFeedback;
extern Pair DiskOnPrimaryAcceleration;
static Pair IndirectTerm;
static real q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static real vt_int[MAX1D], vt_cent[MAX1D];

void ComputeIndirectTerm () {
  IndirectTerm.x = -DiskOnPrimaryAcceleration.x;
  IndirectTerm.y = -DiskOnPrimaryAcceleration.y; 
  if (Indirect_Term == NO) {
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
  }
}
/* Below : work in non-rotating frame */
/* centered on the primary */
void FillForcesArrays (sys)
     PlanetarySystem *sys;
{
  int i, j, l, nr, ns, k, NbPlanets;
  real x, y, angle, distance, distancesmooth;
  real xplanet, yplanet, RRoche,smooth, mplanet;
  real MassTaper, planetmass_taper;
  real PlanetDistance, *Pot, pot, smoothing, *IndPot, *test;
  real InvPlanetDistance3, InvDistance;
  extern boolean MHDLSA, DustFeelDisk, DustFeelPlanets;
  Pot = Potential->Field;
  IndPot = IndPotential->Field;
  test = Test->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  NbPlanets = sys->nb;

  if (MASSTAPER > 1e-2) {
    MassTaper = (PhysicalTime-PhysicalTimeInitial)/(MASSTAPER*2.0*M_PI);
    planetmass_taper = (MassTaper > 1.0 ? 1.0 : .5*(1.0-cos(M_PI*MassTaper)));
  }
  else
    planetmass_taper = 1.0;
  
  /* Indirect term due to the gas acceleration on the star: */
  ComputeIndirectTerm ();

#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) Pot[i] = 0.0;
  for (i = 0; i < (nr+1)*ns; i++) IndPot[i] = 0.0;
  /* -- Gravitational potential from planet on gas -- */
  for (k = 0; k < NbPlanets; k++) {
    xplanet = sys->x[k];
    yplanet = sys->y[k];
    mplanet = sys->mass[k]*planetmass_taper;
    PlanetDistance = sqrt(xplanet*xplanet+yplanet*yplanet);
    InvPlanetDistance3 =  1.0/PlanetDistance/PlanetDistance/PlanetDistance;
    RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
    if (RocheSmoothing)
      smoothing = RRoche*ROCHESMOOTHING;
    else
     smoothing = compute_smoothing (PlanetDistance);
    smooth = smoothing*smoothing;
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,distance,distancesmooth,pot)
    for (i = 0; i < nr; i++) {
      InvDistance = 1.0/Rmed[i];
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	x = Rmed[i]*CosAzimuth[j];
	y = Rmed[i]*SinAzimuth[j];
	distance = (x-xplanet)*(x-xplanet)+(y-yplanet)*(y-yplanet);
	distancesmooth = sqrt(distance+smooth);
	pot = -G*mplanet/distancesmooth; /* Direct term from planet */
	if (Indirect_Term == YES) {
	  pot += G*mplanet*InvPlanetDistance3*(x*xplanet+y*yplanet); /* Indirect term from planet  */
	  if (DustFeelPlanets)
	    IndPot[l] = G*mplanet*InvPlanetDistance3*(x*xplanet+y*yplanet);
	}
	Pot[l] += pot;
      }
    }
  }
  /* -- Gravitational potential from star on gas -- */
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,pot)
  for (i = 0; i < nr; i++) {
    InvDistance = 1.0/Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      x = Rmed[i]*CosAzimuth[j];
      y = Rmed[i]*SinAzimuth[j];
      pot = -G*1.0*InvDistance;  /* Direct term from star */
      /* case where azimuthal extent equals 2pi */
      if ( (!Discard_GasIndirect_term) && (fabs(PMAX-PMIN-2.*M_PI) < 0.01) ) {
	pot -= IndirectTerm.x*x + IndirectTerm.y*y; /* Indirect term from gas on gas */
	if (DustFeelDisk)
	  IndPot[l] -= IndirectTerm.x*x + IndirectTerm.y*y;
	//test[l] = IndPot[l];
      }
      Pot[l] += pot;	
    }
  }
  /* -- Turbulent potential to mimic the MHD turbulence driven by MRI,
     as in Laughlin et al. (2004, LSA04) -- */
  if (MHDLSA)
    ApplyLSAOnPotential();
}

void AdvanceSystemFromDisk (force, Rho, Energy, sys, dt)
     Force *force;
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy;
     real dt;		       
{
  int NbPlanets, k;
  Pair gamma;
  real x, y, r, m, smoothing;
  NbPlanets = sys->nb;
  for (k = 0; k < NbPlanets; k++) {
    if (sys->FeelDisk[k] == YES) {
      m = sys->mass[k];
      x = sys->x[k];
      y = sys->y[k];
      r = sqrt(x*x + y*y);
      if (RocheSmoothing)
	smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
      else
	smoothing = compute_smoothing (r);
      gamma = ComputeAccel (force, Rho, x, y, smoothing, m, sys);
      sys->vx[k] += dt * gamma.x;
      sys->vy[k] += dt * gamma.y;
      /* CB/July 2022: if planet migrates (feels the disc's direct
	 gravity) it should also feel the indirect term arising from
	 the acceleration imprinted by the gas onto the star */
      sys->vx[k] += dt * IndirectTerm.x;
      sys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

void AdvanceSystemRK5 (sys, dt)
     PlanetarySystem *sys;
     real dt;
{
  extern boolean ForcedCircular, ForcedInnerCircular;
  int i, n, myimin;
  boolean *feelothers;
  real dtheta, omega, rdot, x, y, r, v, new_r, vx, vy, theta, denom;
  n = sys->nb;
  //if (ForcedInnerCircular)
  //  myimin = 1;
  //else
  //  myimin = 0; // default case
  if (!ForcedCircular) {
    for (i = 0; i < n; i++) {
      //printf("myimin = %d and i = %d\n",myimin, i);
      q0[i] = sys->x[i];
      q0[i+n] = sys->y[i];
      q0[i+2*n] = sys->vx[i];
      q0[i+3*n] = sys->vy[i];
      PlanetMasses[i] = sys->mass[i];
    }
    feelothers = sys->FeelOthers;
    RungeKunta (q0, dt, PlanetMasses, q1, n, feelothers);
  }
  /* Default case (see below) */
  if (!ForcedInnerCircular) {
    for (i = 1-(PhysicalTime >= RELEASEDATE); i < sys->nb; i++) {
      /* Default case: planets position and velocity updated after 
	 Runge Kutta step */
      if (!ForcedCircular) {
	sys->x[i] = q1[i];
	sys->y[i] = q1[i+n];
	sys->vx[i] = q1[i+2*n];
	sys->vy[i] = q1[i+3*n];
      } else {
	/* Case where planets are held on a fixed circular orbit with 
	   initial angular frequency omega */
	x = sys->x[i];
	y = sys->y[i];
	theta = atan2(y,x);
	vx = sys->vx[i];
	vy = sys->vy[i];
	r = sqrt(x*x + y*y);
	v = sqrt(vx*vx + vy*vy);
	omega = (-y*vx + x*vy)/r/r;
	dtheta = omega*dt;
	sys->x[i]  = r*cos(theta+dtheta);
	sys->y[i]  = r*sin(theta+dtheta);
	sys->vx[i] = -v*sin(theta+dtheta);
	sys->vy[i] =  v*cos(theta+dtheta);
      }
    }
  } else {
    /* New (july 2012): particular case where inner planet held on a fixed 
       circular orbit */
    for (i = 0; i < n; i++) {
      if (i == 0) {  // inner planet (i=0) fixed -> copy-paste of above
	x = sys->x[i];
	y = sys->y[i];
	theta = atan2(y,x);
	vx = sys->vx[i];
	vy = sys->vy[i];
	r = sqrt(x*x + y*y);
	v = sqrt(vx*vx + vy*vy);
	omega = (-y*vx + x*vy)/r/r;
	dtheta = omega*dt;
	sys->x[i]  = r*cos(theta+dtheta);
	sys->y[i]  = r*sin(theta+dtheta);
	sys->vx[i] = -v*sin(theta+dtheta);
	sys->vy[i] =  v*cos(theta+dtheta);
      } else {  // all planets except that indexed with i=0
	sys->x[i] = q1[i];
	sys->y[i] = q1[i+n];
	sys->vx[i] = q1[i+2*n];
	sys->vy[i] = q1[i+3*n];
      }
    }
  }
  /* Case where the innermost planet (with index 0) is drifted
     manually with a prescribed migration rate tuned by RELEASERADIUS 
     and RELEASETIME in .par file */
  if (PhysicalTime < RELEASEDATE) {
    x = sys->x[0];
    y = sys->y[0];
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    rdot = (RELEASERADIUS-r)/(RELEASEDATE-PhysicalTime);
    omega = sqrt((1.+sys->mass[0])/r/r/r);
    new_r = r + rdot*dt;
    denom = r-new_r;
    if (denom != 0.0) {
      dtheta = 2.*dt*r*omega/denom*(sqrt(r/new_r)-1.);
    } else {
      dtheta = omega*dt;
    }
    vx = rdot;
    vy = new_r*sqrt((1.+sys->mass[0])/new_r/new_r/new_r);
    sys->x[0] = new_r*cos(dtheta+theta);
    sys->y[0] = new_r*sin(dtheta+theta);
    sys->vx[0]= vx*cos(dtheta+theta) - vy*sin(dtheta+theta); 
    sys->vy[0]= vx*sin(dtheta+theta) + vy*cos(dtheta+theta); 
  }
}

void SolveOrbits (sys)
     PlanetarySystem *sys;
{
  int i, n;
  real x, y, vx, vy;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    x = sys->x[i];
    y = sys->y[i];
    vx = sys->vx[i];
    vy = sys->vy[i];
    FindOrbitalElements (x, y, vx, vy, 1.0+sys->mass[i], i);
  }
} 

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }
  return lapl;
}

void InitGasDensity (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns, k;
  real *dens, randomnb;
  extern boolean ImposedDensity, AddNoise, AddM1, AddM1Boosted, AddM1toM10;
  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  if (!ImposedDensity) {
    /* Standard case with power-law density profile */
    FillSigma ();
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        dens[l] = SigmaMed[i];
        /* No random noise is added by default to the initial density
          and velocity profiles. If AddNoise set to yes, white noise
          added to the initial density field with arbitrary 1d-3
          relative amplitude by default. This value can be changed with NOISEAMPLITUDE */
        if (AddNoise) {
          randomnb = 2.0*drand48()-1.0;
          dens[l] += NOISEAMPLITUDE*SigmaMed[i]*randomnb;
        }
        if (AddM1) {
          dens[l] += 1e-3*SigmaMed[i]*sin(M_PI*(Rmed[i]-GlobalRmed[0])/(GlobalRmed[GLOBALNRAD-1]-GlobalRmed[0]))*cos(Azimuth[j]);
        }
        if (AddM1Boosted) {
          dens[l] += 1e-1*SigmaMed[i]*sin(M_PI*(Rmed[i]-GlobalRmed[0])/(GlobalRmed[GLOBALNRAD-1]-GlobalRmed[0]))*cos(Azimuth[j]);
        }
        if (AddM1toM10) {
          for (k=0; k<10; k++) {
            dens[l] += 1e-4*SigmaMed[i]*sin(M_PI*(Rmed[i]-GlobalRmed[0])/(GlobalRmed[GLOBALNRAD-1]-GlobalRmed[0]))*cos((k+1.0)*Azimuth[j]);
          }
        }
      }
    }
  } else {
    /* NEW Nov 2010: we impose a fixed density profile */
    InitImposedDensity (Rho);
  }
}


void InitDustDensity (DRho)
     PolarGrid *DRho;
{
  int i, j, l, nr, ns, k;
  real *dens, randomnb;
  extern boolean AddNoise, AddM1, AddM1Boosted, AddM1toM10;
  dens = DRho->Field;
  nr = DRho->Nrad;
  ns = DRho->Nsec;
  FillDSigma ();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      dens[l] = DSigmaMed[i];
      /* No random noise is added by default to the initial density
        and velocity profiles. If AddNoise set to yes, white noise
        added to the initial density field with arbitrary 1d-3
        relative amplitude by default. This value can be changed with NOISEAMPLITUDE */
      if (AddNoise) {
        randomnb = 2.0*drand48()-1.0;
        dens[l] += NOISEAMPLITUDE*DSigmaMed[i]*randomnb;
      }
      if (AddM1) {
	      dens[l] += 1e-3*DSigmaMed[i]*sin(M_PI*(Rmed[i]-GlobalRmed[0])/(GlobalRmed[GLOBALNRAD-1]-GlobalRmed[0]))*cos(Azimuth[j]);
      }
      if (AddM1Boosted) {
	      dens[l] += 1e-1*DSigmaMed[i]*sin(M_PI*(Rmed[i]-GlobalRmed[0])/(GlobalRmed[GLOBALNRAD-1]-GlobalRmed[0]))*cos(Azimuth[j]);
      }
      if (AddM1toM10) {
        for (k=0; k<10; k++) {
          dens[l] += 1e-4*DSigmaMed[i]*sin(M_PI*(Rmed[i]-GlobalRmed[0])/(GlobalRmed[GLOBALNRAD-1]-GlobalRmed[0]))*cos((k+1.0)*Azimuth[j]);
        }
      }
    }
  }
}


void InitImposedDensity (density)
     PolarGrid *density;
{
  int i, ig, j, l, lg, ns, nr;
  real *dens, foo, value;
  FILE *DENSFILE;
  char name_dens[1024];
  real *globaldens;
  sprintf (name_dens, "%saxidens.dat", OUTPUTDIR);
  DENSFILE = fopen (name_dens, "r");
  dens = density->Field;
  nr = density->Nrad;
  ns = density->Nsec;
  if (DENSFILE == NULL) {
    erreur ("ERROR: I could not read the file axidens.dat containing the imposed density profile. Please check and run again\n");
  }
  globaldens = (real*) malloc(sizeof(real)*ns*GLOBALNRAD);
  for (i = 0; i < GLOBALNRAD; i++) {
    fscanf (DENSFILE, "%lf %lf", &foo, &value);
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      globaldens[l] = (real)value;
    }
  }
  for (i = 0; i < nr; i++) {
    ig = (i+IMIN)*ns;
    SigmaMed[i] = globaldens[ig];
    if (i > 0)
      SigmaInf[i] = 0.5*(SigmaMed[i]+SigmaMed[i-1]);
    else
      SigmaInf[i] = SigmaMed[i];
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      lg = (i+IMIN)*ns + j;
      dens[l] = globaldens[lg];
    }
  }
  fclose (DENSFILE);
}

void InitGasEnergy (energ)
     PolarGrid *energ;
{
  int i, j, l, nr, ns;
  real *energy;
  extern boolean TempPresc;
  energy = energ->Field;
  nr = energ->Nrad;
  ns = energ->Nsec;
  FillEnergy ();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energy[l] = EnergyMed[i];
    }
  }
  if (TempPresc)
    FillPrescTime();
}

void InitGasVelocities (Vr, Vt, Rho, DRho)
     PolarGrid *Vr, *Vt, *Rho, *DRho;
{
  extern boolean SGZeroMode;
  extern boolean SelfGravity;
  extern boolean CavityTorque;
  int i, j, l, nr, ns;
  real *vr, *vt, *pres, *cs;
  real *rho, *drho, *St;
  real vk, eta, eps;
  real  r, omega, ri, rii, Hi, myvtheta, vt2;
  real viscosity, t1, t2, r1, r2;
  real num, den;
  real vr_over_cs;
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Vt->Nrad;
  ns  = Vt->Nsec;
  rho = Rho->Field;
  drho = DRho->Field;
  St = Stokes->Field;
  cs = SoundSpeed->Field;
  /* Pressure is already initialized: see initeuler in
     SourceEuler.c */
  pres = Pressure->Field;
  /* ------------------------------------------------------- */
  /* Initialization of azimutal velocity with exact centrifugal
     balance */
  /* ------------------------------------------------------- */
  if ( CentrifugalBalance ) {
    mpi_make1Dprofile (pres, GLOBAL_bufarray);
    /* global axisymmetric pressure field, known by all cpus*/
    for (i = 1; i < GLOBALNRAD; i++) {
      vt_int[i] = ( GLOBAL_bufarray[i] - GLOBAL_bufarray[i-1] ) /	\
	(.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1]) + \
	G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    }
    /* Case with disc self-gravity */
    if ( SelfGravity ) { // Better test with CL rigid!
      if ( !SGZeroMode )
	      mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	      GLOBAL_AxiSGAccr = SG_Accr;
      for (i = 1; i < GLOBALNRAD; i++)
	      vt_int[i] -= ( (Radii[i] - GlobalRmed[i-1])*GLOBAL_AxiSGAccr[i] + \
		       (GlobalRmed[i] - Radii[i])*GLOBAL_AxiSGAccr[i-1] ) / (GlobalRmed[i]-GlobalRmed[i-1]);
    }
    for (i = 1; i < GLOBALNRAD; i++)
      vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
    
    t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
    r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
    t2 = vt_cent[0];
    r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    t1 = t1-r1/(r2-r1)*(t2-t1);
    vt_cent[0] = t1;
    ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[GLOBALNRAD] = vt_cent[GLOBALNRAD-1];
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        vt[l] = vt_cent[i+IMIN];
      }
    }
  } else {
    /* --------- */
    /* Initialization of azimutal velocity without exact centrifugal
       balance (standard procedure): vphi^2 = r/sigma dp/dr + r^2
       Omega_K^2 - ra_{r,sg} where a_{r,sg} denotes the radial
       self-gravitating acceleration */
    /* --------- */
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      omega = sqrt(G*1.0/r/r/r);
      /* default vtheta with density & temperature power-laws */
      myvtheta   = omega*r*sqrt(1.0-pow(ASPECTRATIO,2.0)*		\
				pow(r,2.0*FLARINGINDEX)*		\
				(1.+SIGMASLOPE-2.0*FLARINGINDEX) );
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        vt[l] = myvtheta;
      }
    }
    if (SelfGravity) {
      if ( !SGZeroMode )
	      mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	      GLOBAL_AxiSGAccr = SG_Accr;
      for (i = 0; i <= nr; i++) {
        r = Rmed[i];
        for (j = 0; j < ns; j++) {
          l = i*ns + j;
          vt2 = vt[l]*vt[l];
          /* recall that SG_Accr is centred in radius */
          vt2 -= r*GLOBAL_AxiSGAccr[i+IMIN];
          vt[l] = sqrt(vt2);
        }
      }
    }
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        vt[l] -= OmegaFrame*r;
      }
    }
  }
  /* --------------------------------- */
  /* Initialization of radial velocity */
  /* --------------------------------- */
  if (ViscosityAlpha)
    mpi_make1Dprofile (cs, GLOBAL_SoundSpeed);
  for (i = 0; i < nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    }
    else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    if ((ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0))
      viscosity = FViscosity (r);
    else
      viscosity = 0.0;
    if (CavityTorque)
      vr_over_cs = -( 0.5*(FRACINT+FRACEXT) + 0.5*(FRACINT-FRACEXT)*tanh((CAVITYRADIUS-Rmed[i])/CAVITYWIDTH) );  // negative definite
    for (j = 0; j < ns; j++) {
      l = i*ns+j;
      if (i == nr)
	      vr[l] = 0.0;
      else {
        if (IMPOSEDDISKDRIFT != 0.0)
          vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
        else {
          if (ViscosityAlpha) {
            vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
          } else {
            vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
          }
        }
        if (CavityTorque)
          vr[l] = vr_over_cs*cs[l];  // takes negative values
      }
    }
  }
  /* --------------------------------- */
  /* NEW (04/2020): initialization of dust azimuthal 
     and radial velocities with dust feedback */
  /* --------------------------------- */
  if (DustFluid && DustFeedback) {
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      // Keplerian velocity
      vk = sqrt(G*1.0/r);
      // eta = -1/2 h^2(R) x dlogP / dlogR
      eta = -0.5*pow(ASPECTRATIO,2.0)*pow(r,2.0*FLARINGINDEX)*(-1.-SIGMASLOPE+2.0*FLARINGINDEX);
      // note that GLOBAL_AxiSGAccr has already been computed above
      if (SelfGravity)
	      vk = sqrt(G*1.0/r - r*GLOBAL_AxiSGAccr[i+IMIN]);
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        // eps = dust-to-gas density ratio
        eps = drho[l]/rho[l];
        // gas radial velocity = -eps x dust radial velocity 
        // CUIDADIN! note that centered quantities are used for now...
        vr[l] = 2.0*eps*eta*vk*St[l] / ( St[l]*St[l] + pow(1.0+eps,2.0) );
        // gas azimuthal velocity ~ vk (1-eta) - epsx(dust azimuthal velocity - vk)
        vt[l] = vk*(1.0 - eta) + eps*eta*vk*(1.0+eps) / ( pow(St[l],2.0) + pow(1.0+eps,2.0) );
        vt[l] -= OmegaFrame*r;
      }
    }
  }
  /* NEW (Jan 2015): we keep track of initial radial and azimuthal
     velocities for evanescent boundary function */
  for (i = 0; i < nr; i++) {
    VthetaMed[i] = vt[i*ns]+Rmed[i]*OmegaFrame;
    VradMed[i]   = vr[i*ns];
  }
}

/* CB (04/2020): new initialization that includes dust feedback on gas */
void InitDustVelocities (DVr, DVt, DRho, Rho)
     PolarGrid *DVr, *DVt, *DRho, *Rho;
{
  extern boolean SGZeroMode;
  extern boolean SelfGravity;
  int i, j, l, nr, ns;
  real *drho, *rho, *vr, *vt, *St, *cs, *pres;
  real  r, omega, vk, eta, eps, ri, rii, Hi, myvtheta, vt2, dcs;
  real viscosity, t1, t2, r1, r2;
  real num, den;
  drho= DRho->Field;
  rho = Rho->Field;
  vr  = DVr->Field;
  vt  = DVt->Field;
  nr  = DVt->Nrad;
  ns  = DVt->Nsec;
  St  = Stokes->Field;
  cs  = DSoundSpeed->Field;  
  pres = DPressure->Field;  

  /* ------------------------------------------------------- */
  /* Initialization of azimutal velocity with exact centrifugal
     balance */
  /* ------------------------------------------------------- */
  if ( CentrifugalBalance ) {
    mpi_make1Dprofile (pres, GLOBAL_bufarray);
    /* global axisymmetric pressure field, known by all cpus*/
    for (i = 1; i < GLOBALNRAD; i++) {
      vt_int[i] = ( GLOBAL_bufarray[i] - GLOBAL_bufarray[i-1] ) /	\
	(.5*(DSigma(GlobalRmed[i])+DSigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1]) + \
	G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    }
    /* Case with disc self-gravity */
    if ( SelfGravity ) { // Better test with CL rigid!
      if ( !SGZeroMode )
	      mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	      GLOBAL_AxiSGAccr = SG_Accr;
      for (i = 1; i < GLOBALNRAD; i++)
	      vt_int[i] -= ( (Radii[i] - GlobalRmed[i-1])*GLOBAL_AxiSGAccr[i] + \
		       (GlobalRmed[i] - Radii[i])*GLOBAL_AxiSGAccr[i-1] ) / (GlobalRmed[i]-GlobalRmed[i-1]);
    }
    for (i = 1; i < GLOBALNRAD; i++)
      vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
    
    t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
    r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
    t2 = vt_cent[0];
    r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    t1 = t1-r1/(r2-r1)*(t2-t1);
    vt_cent[0] = t1;
    ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[GLOBALNRAD] = vt_cent[GLOBALNRAD-1];
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        vt[l] = vt_cent[i+IMIN];
      }
    }
  } else {
  /* --------- */
  /* Initialization of azimutal velocity without exact centrifugal
     balance (standard procedure): vphi^2 = r/sigma dp/dr + r^2
     Omega_K^2 - ra_{r,sg} where a_{r,sg} denotes the radial
     self-gravitating acceleration */
  /* --------- */
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      omega = sqrt(G*1.0/r/r/r);
      /* default vtheta with density & temperature power-laws */
      myvtheta   = omega*r*sqrt(1.0-pow(DASPECTRATIO,2.0)*		\
				pow(r,2.0*DFLARINGINDEX)*		\
				(1.+SIGMASLOPE-2.0*DFLARINGINDEX) );
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        vt[l] = myvtheta;
      }
    }
    if (SelfGravity) {
      if ( !SGZeroMode )
	      mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	      GLOBAL_AxiSGAccr = SG_Accr;
      for (i = 0; i <= nr; i++) {
        r = Rmed[i];
        for (j = 0; j < ns; j++) {
          l = i*ns + j;
          vt2 = vt[l]*vt[l];
          /* recall that SG_Accr is centred in radius */
          vt2 -= r*GLOBAL_AxiSGAccr[i+IMIN];
          vt[l] = sqrt(vt2);
        }
      }
    }
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        vt[l] -= OmegaFrame*r;
      }
    }
  }
  /* --------------------------------- */
  /* Initialization of radial velocity */
  /* --------------------------------- */
  if (DViscosityAlpha)
    mpi_make1Dprofile (cs, GLOBAL_DustSoundSpeed);
  for (i = 0; i < nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    }
    else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    if ((DALPHAVISCOSITY != 0.0) || (DVISCOSITY != 0.0))
      viscosity = DFViscosity (r);
    else
      viscosity = 0.0;
    for (j = 0; j < ns; j++) {
      l = i*ns+j;
      if (i == nr)
	      vr[l] = 0.0;
      else {
	      vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
        if (DViscosityAlpha) {
          vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*DFLARINGINDEX+1.0);
        } else {
          vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
        }
      }
    }
  }
  /* --------------------------------- */
  /* NEW (04/2020): initialization of dust azimuthal 
     and radial velocities with dust feedback */
  /* --------------------------------- */
  if (DustFeedback) {
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      vk = sqrt(G*1.0/r);
      // eta = -1/2 h^2(R) x dlogP / dlogR
      eta = -0.5*pow(ASPECTRATIO,2.0)*pow(r,2.0*FLARINGINDEX)*(-1.-SIGMASLOPE+2.0*FLARINGINDEX);
      // note that GLOBAL_AxiSGAccr has already been computed above...
      if (SelfGravity)
	      vk = sqrt(G*1.0/r - r*GLOBAL_AxiSGAccr[i+IMIN]);
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        // eps = dust-to-gas density ratio
        eps = drho[l]/rho[l];
        // radial dust velocity
        vr[l] = -2.0*eta*vk / ( St[l] + pow(1.0+eps,2.0)/St[l] );
        // azimuthal dust velocity
        vt[l] = vk - eta*vk*(1.0+eps) / ( pow(St[l],2.0) + pow(1.0+eps,2.0) );
        vt[l] -= OmegaFrame*r;
      }
    }
  } else {
    for (i = 0; i < nr; i++) {
      r = Rmed[i];
      vk = sqrt(G*1.0/r);
      // eta = -1/2 h^2(R) x dlogP / dlogR
      eta = -0.5*pow(ASPECTRATIO,2.0)*pow(r,2.0*FLARINGINDEX)*(-1.-SIGMASLOPE+2.0*FLARINGINDEX);
      // note that GLOBAL_AxiSGAccr has already been computed above...
      if (SelfGravity)
        vk = sqrt(G*1.0/r - r*GLOBAL_AxiSGAccr[i+IMIN]);
      for (j = 0; j < ns; j++) {
        l = i*ns + j;
        // radial dust velocity
        vr[l] = -2.0*eta*vk / ( St[l] + pow(St[l],-1.0) );
        // azimuthal dust velocity
        vt[l] = vk - eta*vk / ( pow(St[l],2.0) + 1.0);
        vt[l] -= OmegaFrame*r;
      }
    }
  }

  /* NEW (Jan 2015): we keep track of initial radial and azimuthal
     velocities for evanescent boundary function */
  for (i = 0; i < nr; i++) {
    DVthetaMed[i] = vt[i*ns]+Rmed[i]*OmegaFrame;
    DVradMed[i]   = vr[i*ns];
  }
}


/* CB (04/2020): initialise dust Stokes number assuming the relative
   velocity between gas and dust is much smaller than the sound speed,
   which is generally expected as long as the dust/gas density ratio
   takes reasonable values (ie not much larger than unity) */
void InitDustStokesNumber (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  real Cdrag, k_D, f_D, Kn, lambda, densg_cgs;
  real dust_internal_density, dust_radius;
  real *dens, *St;

  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  St   = Stokes->Field;

  /* Convert particle's mean density from g/cm^3 to code units */
  dust_internal_density = RHOPART*1000.0*pow(unit_length, 3.)/unit_mass;
  /* Convert particle's radius from meters to code units */
  dust_radius = Sizepart/unit_length;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* Molecular mean-free path ~ m_H2 H / Sigma_gas d^2 with m_H2
	 mass of H2, H pressure scale height and d ~ typical
	 diameter of a particle */
      densg_cgs = 0.1*dens[l]*unit_mass*pow(unit_length,-2.0);  // in g cm^-2
      lambda = 3.34e-8 * pow(densg_cgs,-1.) * AspectRatio(Rmed[i]) * Rmed[i]; // in code units
      /* Knudsen number, note that dust_radius is already in code units */
      Kn = 0.5*lambda/dust_radius;
      /* f_D coefficient to link the subsonic and the supersonic regimes */
      f_D = 1.0;
      /* k_D coefficient that describes the Stokes regime */
      k_D = 1.0;
      /* Final expression for drag coefficient Cdrag as in Paardekooper 07 */
      Cdrag = pow(3.0*Kn+1.0,2.0) * pow(9.0*Kn*Kn*f_D + 3.0*Kn*k_D,-1.0);
      /* Dimensionless stopping time at the particle = stokes number */
      St[l] =  0.5*M_PI*Cdrag*dust_internal_density*dust_radius/dens[l];
    }
  }
  
}
