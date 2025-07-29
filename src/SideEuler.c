/** \file SideEuler.c

    Total mass and angular momentum monitoring, and boundary conditions.
    In addition, this file contains a few low-level functions that
    manipulate PolarGrid 's or initialize the forces evaluation.

*/

#include "mp.h"

extern boolean OpenInner, OpenInnerDust, KNOpen, NonReflecting, OuterSourceMass, Evanescent;
extern boolean SelfGravity, SGZeroMode, EnergyEquation, MixedBC, AccBoundary;
extern Pair DiskOnPrimaryAcceleration;
real Hp0, Hg0, Ht0;

real GasTotalMass (array)
     PolarGrid *array;
{
  int i, j, ns;
  real *density, total = 0.0, fulltotal=0.0;
  ns = array->Nsec;
  density = array->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      total += Surf[i]*density[j+i*ns];
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}

real GasMomentum (Density, Vtheta)
     PolarGrid *Density, *Vtheta;
{
  int i,j,l,ns;
  real vt_cent;
  real *density, *vtheta, total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vtheta = Vtheta->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      /* centered-in-cell azimuthal velocity */
      if (j < ns-1)
	vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
	vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      total += Surf[i]*density[l]*Rmed[i]*vt_cent;
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}

real GasTotalEnergy (Density, Vrad, Vtheta, Energy)
     PolarGrid *Density, *Vrad, *Vtheta, *Energy;
{
  int i, j, l, ns;
  real *density, *vrad, *vtheta, *energy, *pot;
  real vr_cent, vt_cent;
  real total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  energy = Energy->Field;
  pot = Potential->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* centered-in-cell radial velocity */
      vr_cent = (Rmed[i]-Rinf[i])*vrad[l+ns] + (Rsup[i]-Rmed[i])*vrad[l];
      vr_cent /= (Rsup[i]-Rinf[i]);
      /* centered-in-cell azimuthal velocity */
      if (j < ns-1)
	vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
	vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      total += 0.5*Surf[i]*density[l]*(vr_cent*vr_cent + vt_cent*vt_cent) + \
	Surf[i]*energy[l] -						\
	Surf[i]*density[l]*pot[l];
      /* Gas total energy is the sum of its kinematic energy, internal energy */
      /* and gravitational potential energy, including self-gravity */
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}


void CheckMomentumConservation (Density, Vtheta, sys)
     PolarGrid *Density, *Vtheta;
     PlanetarySystem *sys;
{
  FILE *fichmom;
  char name[256];
  int k;
  real totalmomentum, plmom;
  real xplanet, yplanet, vxplanet, vyplanet;
  real rpl, thetapl, vazimpl, masspl;
  real gasmom, planetsmom;
  gasmom = GasMomentum (Density, Vtheta);
  planetsmom = 0.;
  
  for ( k = 0; k < sys->nb; k++ ) {
    xplanet     = sys->x[k];
    yplanet     = sys->y[k];
    rpl         = sqrt( xplanet*xplanet + yplanet*yplanet );
    thetapl     = atan2 (yplanet, xplanet);
    vxplanet    = sys->vx[k];
    vyplanet    = sys->vy[k];
    vazimpl     = -vxplanet*sin(thetapl) + vyplanet*cos(thetapl);
    masspl      = sys->mass[k];
    plmom       = masspl*rpl*vazimpl;
    planetsmom += plmom;
  }
  totalmomentum = gasmom + planetsmom;
  if ( PhysicalTime < 1e-10 ) {
    Hp0 = plmom;
    Hg0 = gasmom;
    Ht0 = totalmomentum;
    printf("time = %lg, Hp0 = %lg, Hg0 = %lg et Ht0 = %lg\n", PhysicalTime, Hp0, Hg0, Ht0);
  }
  if (!CPU_Master) return;
  sprintf (name, "%s%s.dat", OUTPUTDIR, "Momentum");
  fichmom = fopen(name, "a");
  if (fichmom == NULL) {
    fprintf (stderr, "Can't write 'Momentum.dat' file. Aborting.\n");
    prs_exit (1);
  }
  plmom = fabs (plmom - Hp0);
  gasmom = fabs (gasmom - Hg0);
  totalmomentum = fabs (totalmomentum - Ht0);
  fprintf (fichmom, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", PhysicalTime, plmom, gasmom, totalmomentum, totalmomentum / Ht0);
  fclose (fichmom);
}


void DivisePolarGrid (Num, Denom, Res)
     PolarGrid *Num, *Denom, *Res;
{
  int i,j,l,nr,ns;
  real *num, *denom, *res;
  num = Num->Field;
  denom=Denom->Field;
  res = Res->Field;
  ns = Res->Nrad;
  nr = Res->Nsec;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      res[l] = num[l]/(denom[l]+1e-20);
    }
  }
}

void InitComputeAccel ()
{
  int i, j, l, nr, ns;
  real *abs, *ord;
  CellAbscissa = CreatePolarGrid (NRAD,NSEC,"abscissa");
  CellOrdinate = CreatePolarGrid (NRAD,NSEC,"ordinate");
  nr = CellAbscissa->Nrad;
  ns = CellAbscissa->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      abs[l] = Rmed[i] * cos(Azimuth[j]);
      ord[l] = Rmed[i] * sin(Azimuth[j]);
    }
  }
}
  
Pair ComputeAccel (force, Rho, x, y, rsmoothing, mass, psys)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
     PlanetarySystem *psys;
{
  Pair acceleration;
  ComputeForce (force, Rho, x, y, rsmoothing, mass, psys);
  if (ExcludeHill) {
    acceleration.x = force->fx_ex_inner+force->fx_ex_outer;
    acceleration.y = force->fy_ex_inner+force->fy_ex_outer;
  } else {
    acceleration.x = force->fx_inner+force->fx_outer;
    acceleration.y = force->fy_inner+force->fy_outer;
  }
  return acceleration;
}

void OpenBoundary (Vrad, Vtheta, Rho, Energy)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
{
  int i,j,l,ns,nr;
  extern boolean DontApplySubKeplerian, RadiativeDiffusion, CavityTorque;
  real *rho, *vr, *vt, *energy, *cs, visc, vri, *kcoeff;
  real R0, R1, R2, Rsup0, Rsup1, Ksup0, Ksup1, T0, T1, T2;
  real Rnm3, Rnm2, Rnm1, Rsupnm3, Rsupnm2, Ksupnm3, Ksupnm2, Tnm3, Tnm2, Tnm1;
  cs = SoundSpeed->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  vt  = Vtheta->Field;
  energy = Energy->Field;
  kcoeff = RadiativeKCoeff->Field;
  /* -------------------------------- */
  /* Inner Boundary Condition         */
  /* -------------------------------- */
  if (CPU_Rank == 0) {
    i = 2;  // NEW July 2023
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (KNOpen) {
        /* Kley and Nelson (2008) prescription */
        rho[l-ns] = rho[l]*Rmed[1]/Rmed[0];  /* copy first ring into ghost ring */
        if (EnergyEquation)
          energy[l-ns] = energy[l]; // zero gradient for thermal energy
        if ((vr[l+ns] > 0.0))
          vr[l] = 0.0; /* we just allow outflow [inwards] */
        else {
          vri = -3.*FViscosity(Rmed[i])/Rmed[i]*3.*(-SIGMASLOPE+1.+2.*FLARINGINDEX); 
          if(fabs(vr[l+ns])>fabs(vri))
            vr[l]=vri;
          else
            vr[l]=vr[l+ns];
        }
      } else {
        /* Standard outflow prescription */
        /*
        if (ViscosityAlpha || (VISCOSITY != 0.0) ) // set nuxSigma uniform
          rho[l-ns] = rho[l] * FViscosity(Rmed[i]) / FViscosity(Rmed[i-1]);
          vr[l] = -1.5*FViscosity(Rinf[i])/Rinf[i];
        else
          rho[l-ns] = rho[l] ;      // zero gradient for surface density
          vr[l] = 0.0;
        */
        rho[l-ns]   = rho[l];      // zero gradient for surface density   // cuidadin
        rho[l-2*ns] = rho[l];      // zero gradient for surface density   // cuidadin
        if (EnergyEquation) {
          energy[l-ns]   = energy[l]; // zero gradient for thermal energy
          energy[l-2*ns] = energy[l]; // zero gradient for thermal energy
        }
        if (DontApplySubKeplerian) {
          vt[l-ns]   = vt[l] * sqrt(Rmed[i]/Rmed[i-1]);  // Keplerian extrapolation
          vt[l-2*ns] = vt[l] * sqrt(Rmed[i]/Rmed[i-2]);  // Keplerian extrapolation
        }
        // we just allow inflow
        if ((vr[l+ns] >= 0.0) || (rho[l] < SigmaMed[0])) {
          vr[l]    = 0.0;
          vr[l-ns] = 0.0;
        }
        else {
          vr[l]    = vr[l+ns];
          vr[l-ns] = vr[l+ns];
        }
        /*
        else { // openinner
          if ((vr[l+ns] >= 0.0) || (rho[l] < SigmaMed[0]))
            vr[l] = 0.0;            // vr set to zero when directed outward
          else
            vr[l] = vr[l+ns];       // vr extrapolated otherwise 
        }
        */
      }
      if (CavityTorque) {
        /* June 2022 */
        rho[l-ns] = SigmaMed[0];
        vt[l-ns]  = VthetaMed[0]-Rmed[0]*OmegaFrame;;
        vr[l]     = VradMed[0];
      }
    }
  }
  /* -------------------------------- */
  /* Outer boundary condition         */
  /* -------------------------------- */
  //if ( (CPU_Rank == CPU_Highest) && (!MixedBC) ) { // cuidadin
  if ( CPU_Rank == CPU_Highest ) {
    i = nr-1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /*
      if (ViscosityAlpha || (VISCOSITY != 0.0) )
	      rho[l] = rho[l-ns] * FViscosity(Rmed[i-1]) / FViscosity(Rmed[i]);
      else
	      rho[l] = rho[l-ns];	  // zero gradient for surface density
      */
      rho[l] = rho[l-ns];
      if (EnergyEquation)
	      energy[l] = energy[l-ns]; // zero gradient for thermal energy
      if (DontApplySubKeplerian)
	      vt[l] = vt[l-ns] * pow(Rmed[i]/Rmed[i-1],-0.5);
      if (KNOpen) {
        if ((vr[l-ns] < 0.0))
          vr[l]  = 0.0; // we just allow outflow
        else
          vr[l] = vr[l-ns];
      }
      else  {
        /*
        if (ViscosityAlpha || (VISCOSITY != 0.0) )
          vr[l] = -1.5*FViscosity(Rinf[i])/Rinf[i];
        else
          vr[l] = 0.0;
        */
        //vr[l] = vr[l-ns]; // CUIDADIN!!
        // Below BC leads to oscillating behavior in vr... 
        if ((vr[l-ns] < 0.0) || (rho[l] < SigmaMed[nr-2]))
          vr[l] = 0.0;            // vr set to zero when directed inward
        else
          vr[l] = vr[l-ns];       // vr extrapolated otherwise 
      }
      if (CavityTorque) {
        /* June 2022  */
        rho[l] = SigmaMed[i];
        vt[l]  = VthetaMed[i]-Rmed[i]*OmegaFrame;
        vr[l]  = VradMed[i];
      }
    }
  }
}


void OpenBoundaryd (DVrad, DVtheta, DRho)
     PolarGrid *DVrad, *DVtheta, *DRho;
{
  int i,j,l,ns,nr;
  extern boolean DontApplySubKeplerian;
  real *rho, *vr, *vt, visc, vri;
  ns = DRho->Nsec;
  nr = DRho->Nrad;
  rho = DRho->Field;
  vr  = DVrad->Field;
  vt  = DVtheta->Field;
  /* -------------------------------- */
  /* Inner Boundary Condition         */
  /* -------------------------------- */
  if (CPU_Rank == 0) {
    i = 2;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* Standard outflow prescription */
      rho[l-ns] = rho[l] ;      // zero gradient for surface density
      rho[l-2*ns] = rho[l];      // zero gradient for surface density   // cuidadin
      /* NEW (Sept 28 2011): if subkeplerian BC is not applied, then
        do power-law extrapolation */
      if (DontApplySubKeplerian) {
        vt[l-ns] = vt[l] * pow(Rmed[i-1]/Rmed[i],-0.5);
        vt[l-2*ns] = vt[l] * pow(Rmed[i-2]/Rmed[i],-0.5);
      }
      /*
      if ((vr[l+ns] >= 0.0) || (rho[l] < DSigmaMed[0]))
        vr[l] = 0.0;            // vr set to zero when directed outward
      else
        vr[l] = vr[l+ns];       // vr extrapolated otherwise 
      */
      vr[l] = vr[l+ns]; // CUIDADIN!
      vr[l-ns] = vr[l+ns]; // CUIDADIN!
    }
  }
  /* -------------------------------- */
  /* Outer boundary condition         */
  /* -------------------------------- */
  if ( (CPU_Rank == CPU_Highest) && (!MixedBC) ) {
    i = nr-1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rho[l] = rho[l-ns];	// zero gradient for surface density
      /* NEW (Sept 28 2011): if subkeplerian BC is not applied, then
	   do power-law extrapolation */
      if (DontApplySubKeplerian)
	      vt[l] = vt[l-ns] * pow(Rmed[i]/Rmed[i-1],-0.5);
      /*
      if ((vr[l-ns] < 0.0) || (rho[l] < DSigmaMed[nr-2]))
        vr[l] = 0.0;            // vr set to zero when directed inward
      else
        vr[l] = vr[l-ns];       // vr extrapolated otherwise 
      */
	    vr[l] = vr[l-ns]; // CUIDADIN!
    }
  }
}


/* New BC developed by Sareh Ataiee for a disk with imposed,
   time-varying Mdot  */
void AccretingBoundary (Vrad, Vtheta, Rho, Energy)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
{
  int i,j,l,ns,nr;
  real *dens, *vr, *vt, *energy;
  real mdotinit, mdotfinal, mdot;
  real r, vrm, vrp;
  real omega, vr_med, axi[GLOBALNRAD];
  extern boolean DecInner, Restart;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  dens = Rho->Field;
  vr  = Vrad->Field;
  vt  = Vtheta->Field;
  energy = Energy->Field;
  /* We'll need the global axisymmetric density profile */
  mpi_make1Dprofile (dens, axi);
  /* -------------------------------- */
  /* Inner boundary condition         */
  /* -------------------------------- */
  if (CPU_Rank == 0) {
    i = 1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energy[l-ns] = energy[l];
      /* Option to decrease surface density in first ring */
      if (DecInner){
	if (PhysicalTime != 0.0)
	  dens[l] *= 0.999;
	if (vr[l+ns] > 0)
	  vr[l] = 0;
	else
	  vr[l] = vr[l+ns];
      } else {
	vr[l] = -1.5*FViscosity(Rinf[i])/Rinf[i];
	// Line below not strictly necessary (testing purposes mainly)
	vr[l-ns] = -1.5*FViscosity(Rinf[i-1])/Rinf[i-1];
	/* Zero gradient in Mdot ~ nu Sigma */
	dens[l-ns] = dens[l] * FViscosity(Rinf[i]) / FViscosity(Rinf[i-1]);
	vt[l-ns] = VthetaMed[i-1]-Rmed[i-1]*OmegaFrame;
      }
    }
  }
  /* -------------------------------- */
  /* Outer boundary condition         */
  /* -------------------------------- */
  if (CPU_Rank == CPU_Highest) {
    i = nr-1;
    r = GlobalRmed[GLOBALNRAD-1];
    vrm = -1.5*FViscosity(Radii[GLOBALNRAD-1])/Radii[GLOBALNRAD-1];
    vrp = -1.5*FViscosity(Radii[GLOBALNRAD])/Radii[GLOBALNRAD];
    /* We calculate Mdot so that it's consistent with the initial 
       gas surface density and radial velocity at the outer edge: 
       Mdot = -2pi r_out Sigma(r_out) vr(r_out) */
    mdotinit = fabs(2.0*PI*Sigma(r)*r*0.5*(vrp+vrm));
    mdotfinal = mdotinit; // CUIDADIN
    /* We slowly decrease Mdot by one order of magnitude (default pace)
       over a user-prescribed timescale, called MDOTTIME */
    if ((PhysicalTime/2./PI) <= MDOTTIME){
      mdot = mdotinit + (mdotfinal - mdotinit) * (PhysicalTime/2./PI) / MDOTTIME;
      mdot *= -1.0;
    } else {
      mdot = -mdotfinal;
    }
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;      
      energy[l] = energy[l-ns];
      vr[l] = -1.5*FViscosity(Rinf[i])/Rinf[i];
      // Line below not strictly necessary (testing purposes mainly)
      vr[l+ns] = -1.5*FViscosity(Rinf[i+1])/Rinf[i+1];
      vr_med = 0.5*(vr[l]+vr[l+ns]);
      dens[l] = mdot/2./PI/Rmed[i]/vr_med;
      vt[l] = VthetaMed[i]-Rmed[i]*OmegaFrame;
      //printf ("vr_med = %lg, mdot = %lg, dens = %lg, vt = %lg\n",vr_med,mdot,dens[l],vt[l]);
    }
    if (axi[i+IMIN-2] <= floordens){
      for (j = 0; j < ns; j++){
	l = i*ns+j;
	dens[l] = floordens;
	dens[l-ns] = floordens;
      }
    } else if (axi[i+IMIN-1] <= floordens){
      for (j = 0; j < ns; j++){
	l = i*ns+j;
	dens[l] = floordens;
      }
    }
  }
}


void NonReflectingBoundary (Vrad, Rho, Energy, Vtheta, sys)
     PolarGrid *Vrad, *Rho, *Energy, *Vtheta;
     PlanetarySystem *sys;
{
  int i, j, k, l, ns, nr, jp, lp, i_angle;
  real *rho, *vr, *vt, *cs, *energy;
  real dangle, mean;
  real cs0, cs1, vr_med, csnrm1, csnrm2;
  real dens0, dens1, densnrm2, densnrm1;
  real csinit0, csinit1;
  real xp, yp, rp, r_in, omega_in, r_out, omega_out, buf;
  cs = SoundSpeed->Field;
  energy = Energy->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  vt = Vtheta->Field;
  /* ======================================= */
  /* INNER NON-REFLECTING BOUNDARY CONDITION */
  /* ======================================= */
  if (CPU_Rank == 0) {
    i=1;
    /* Find radius and angular frequency of innermost planet */
    buf = Radii[GLOBALNRAD];
    for (k=0; k<sys->nb; k++) {
      xp = sys->x[k];
      yp = sys->y[k];
      rp = sqrt( xp*xp + yp*yp );
      if ( (rp < buf) && (rp >= Radii[0]) ) {
	r_in = rp;
	omega_in = pow(rp,-1.5);
      }
      buf = r_in;
    }
    /* ------------------------ */
    /* STANDARD BAROTROPIC CASE */
    /* ------------------------ */
    if (!EnergyEquation) {
      /* (a) sound speed in first two rings */
      cs0 = 0.0;
      cs1 = 0.0;
      dens0 = 0.0;
      dens1 = 0.0;
      for (j = 0; j < ns; j++) {
	cs0 += cs[j];
	cs1 += cs[ns+j];
	dens0 += rho[j];
	dens1 += rho[ns+j];
      }
      cs0 /= (real)ns;
      cs1 /= (real)ns;
      dens0 /= (real)ns;
      dens1 /= (real)ns;
      /* (b) pitch angle calculation */
      dangle = (pow(Rinf[1],-1.5)-omega_in)/(.5*(cs0+cs1));
      dangle *= (Rmed[1]-Rmed[0]);
      i_angle = (int)(dangle/(PMAX-PMIN)*(real)NSEC+.5);
      /* (c) boundary on vrad and dens fields */
#pragma omp parallel for private(l,jp,lp,vr_med)
      for (j = 0; j < ns; j++) {
	l = j+i*ns;   // recall i=1
	jp = j+i_angle;
	if (jp >= ns) jp -= ns;
	if (jp < 0) jp += ns;
	lp = jp;
	if (AccBoundary == NO)
	  rho[lp] = rho[l];       /* copy first ring into ghost ring */
	else
	  rho[lp] = dens0 + (rho[l]-dens1);
	vr_med = -cs1*(rho[l]-dens1)/dens1;
	if (AccBoundary == YES) {
	  vr[l] += (2.*vr_med-vr[l+ns]);
	}
	else
	  vr[l] = 2.*vr_med-vr[l+ns];
      }
      if (AccBoundary == NO) {
	/* (d) density adjustment */
	mean = 0.0;
	for (j = 0; j < ns; j++) {
	  mean += rho[j];
	}
	mean /= (real)ns;
	for (j = 0; j < ns; j++) {
	  rho[j] += (SigmaMed[0]-mean);
	}
      }
    }
    /* ------------------------ */
    /*      ADIABATIC CASE      */
    /* ------------------------ */
    if (EnergyEquation) {
      /* (a) sound speed in first two rings */
      cs0 = 0.0;
      cs1 = 0.0;
      for (j = 0; j < ns; j++) {
	cs0 += cs[j];
	cs1 += cs[ns+j];
      }
      cs0 /= (real)ns;
      cs1 /= (real)ns;
      csinit0 = sqrt(ADIABATICINDEX)*ASPECTRATIO*pow(Rmed[0],-0.5+FLARINGINDEX);
      csinit1 = sqrt(ADIABATICINDEX)*ASPECTRATIO*pow(Rmed[1],-0.5+FLARINGINDEX);
      dangle = (pow(Rinf[1],-1.5)-omega_in)/(.5*(csinit0+csinit1));
      dangle *= (Rmed[1]-Rmed[0]);
      i_angle = (int)(dangle/(PMAX-PMIN)*(real)NSEC+.5);
      /* (a) shift on energy and vrad */
#pragma omp parallel for private(l,jp,lp,vr_med)
      for (j = 0; j < ns; j++) {
	/* The expression below should be refined as we need to know the
	   orbital frequency of the nearest planet */
	l = j+i*ns;   // recall i=1
	jp = j+i_angle;
	if (jp >= ns) jp -= ns;
	if (jp < 0) jp += ns;
	lp = jp;
	energy[lp] = energy[l];
	vr_med = -csinit1*(rho[l]-SigmaMed[1])/SigmaMed[1];
	vr[l] = 2.*vr_med-vr[l+ns];
      }
      /* (b) shift on acoustic part of density */
#pragma omp parallel for private(l,jp,lp,vr_med)
      for (j = 0; j < ns; j++) {
	/* The expression below should be refined as we need to know the
	   orbital frequency of the nearest planet */
	rho[j] = SigmaMed[0] + (ADIABATICINDEX-1.0)*(energy[j]-EnergyMed[0])/csinit0/csinit0;
      }
      /* (c) density and energy adjustments */
      mean = 0.0;
      for (j = 0; j < ns; j++) {
	mean += rho[j];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
	rho[j] += (SigmaMed[0]-mean);
      }
      mean = 0.0;
      for (j = 0; j < ns; j++) {
	mean += energy[j];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
	energy[j] += (EnergyMed[0]-mean);
      }
    }
  }
  /* ======================================= */
  /* OUTER NON-REFLECTING BOUNDARY CONDITION */
  /* ======================================= */
  if (CPU_Rank == CPU_Highest) {
    /* Find radius and angular frequency of outermost planet */
    buf = Radii[0];
    for (k=0; k<sys->nb; k++) {
      xp = sys->x[k];
      yp = sys->y[k];
      rp = sqrt( xp*xp + yp*yp );
      if ( (rp > buf) && (rp <= Radii[GLOBALNRAD]) ) {
	r_out = rp;
	omega_out = pow(rp,-1.5);
      }
      buf = r_out;
    }
    csnrm2 = 0.0;
    csnrm1 = 0.0;
    densnrm1 = 0.0;
    densnrm2 = 0.0;
    for (j=0; j<ns; j++) {
      csnrm2 += cs[(nr-2)*ns+j];
      csnrm1 += cs[(nr-1)*ns+j];
      densnrm1 += rho[(nr-1)*ns+j];
      densnrm2 += rho[(nr-2)*ns+j];
    }
    csnrm1 /= (real)ns;
    csnrm2 /= (real)ns;
    densnrm1 /= (real)ns;
    densnrm2 /= (real)ns;
    i = nr-1;		 /* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[nr-2],-1.5)-omega_out)/(.5*(csnrm1+csnrm2));
    dangle *= (Rmed[nr-1]-Rmed[nr-2]);
    i_angle = (int)(dangle/(PMAX-PMIN)*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j-i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp+(i-1)*ns;
      if (AccBoundary == NO) {
	rho[l] = rho[lp];		/* copy first ring into ghost ring */
	energy[l] = energy[lp];	        /* copy first ring into ghost ring */
      } else {
	rho[l] = densnrm1 + (rho[lp]-densnrm2);
	energy[l] = energy[lp];	        /* copy first ring into ghost ring */
      }
      if (!EnergyEquation)
	vr_med = csnrm1*(rho[l-ns]-densnrm2)/densnrm2;
      else
	vr_med = cs[l]*(rho[l-ns]-densnrm2)/densnrm2;
      if (AccBoundary == YES) {
	vr[l] += (2.*vr_med-vr[l-ns]);
      }
      else
	vr[l] = 2.*vr_med-vr[l-ns];
    }
    /* density and energy adjustments */
    if (AccBoundary == NO) {
      mean = 0.0;
      for (j = 0; j < ns; j++) {
	mean += rho[j+ns*(nr-1)];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
	rho[j+(nr-1)*ns] += SigmaMed[nr-1]-mean;
      }
      mean = 0.0;
      for (j = 0; j < ns; j++) {
	mean += energy[j+ns*(nr-1)];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
	energy[j+(nr-1)*ns] += EnergyMed[nr-1]-mean;
      }
    }
  }
}

void EvanescentBoundary (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, step)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     PolarGrid *DVrad, *DVtheta, *DRho;
     real step;
{
  int i, j, l, nr, ns;
  int myi1D;
  real *vrad, *vtheta, *dens, *energ;
  real *dvrad, *dvtheta, *ddens;
  real vrad0, vtheta0, dens0, energ0;
  real dvrad0, dvtheta0, ddens0;
  real damping, Tin, Tout, lambda;
  extern boolean DampToIni, DampToAxi, DampToViscous, DustFluid;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  if (DustFluid && (OpenInnerDust == NO)) {
    dvrad = DVrad->Field;
    dvtheta = DVtheta->Field;
    ddens = DRho->Field;
  }
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  
  lambda = 0.0;
  for (i = 0; i < nr; i++) {
    if ( (Rmed[i] < WKZRMIN) || (Rmed[i] > WKZRMAX) ) {
      /* Damping operates only inside the wave killing zones */
      if (Rmed[i] < WKZRMIN) {
	damping = (Rmed[i]-WKZRMIN)/(GlobalRmed[0]-WKZRMIN);
	Tin = 0.3*pow(Rmed[i],1.5);  // 0.3 Omega_K^-1
	lambda = damping*damping*step/Tin;
      }
      if (Rmed[i] > WKZRMAX) {
	damping = (Rmed[i]-WKZRMAX)/(GlobalRmed[GLOBALNRAD-1]-WKZRMAX);
	Tout = 0.3*pow(Rmed[i],1.5);  // 0.3 Omega_K^-1
	lambda = damping*damping*step/Tout;
      }
      /* Damping wrt initial profiles (requires DampToIni set to yes
	 in .par file) */
      if (DampToIni) {
	vtheta0 = VthetaMed[i]-Rmed[i]*OmegaFrame;
	vrad0 = VradMed[i];
	dens0 = SigmaMed[i];
	energ0 = EnergyMed[i];
      }
      if (DustFluid && DampToIni && (OpenInnerDust == NO)) {
	dvtheta0 = DVthetaMed[i]-Rmed[i]*OmegaFrame;
	dvrad0 = DVradMed[i];
	ddens0 = DSigmaMed[i];
      }
      /* damping wrt instantaneous axisymmetric fields (default
	 case) */
      if (DampToAxi) {
	vrad0   = 0.0;
	vtheta0 = 0.0;
	dens0   = 0.0;
	energ0  = 0.0;
	for (j = 0; j < ns; j++) {
	  l = i*ns + j;
	  vrad0   += vrad[l];
	  vtheta0 += vtheta[l];
	  dens0   += dens[l];
	  energ0  += energ[l];
	}
	vrad0   /= (real)ns;
	vtheta0 /= (real)ns;
	dens0   /= (real)ns;
	energ0  /= (real)ns;
      }
      /* Damping wrt viscous evolving 1D profiles (requires
	 DampToViscous set to yes in .par file) */
      if (DampToViscous) {
	vtheta0 = 0.0;
	energ0  = 0.0;
	for (j = 0; j < ns; j++) {
	  l = i*ns + j;
	  vtheta0 += vtheta[l];
	  energ0  += energ[l];
	}
	vtheta0 /= (real)ns;
	energ0  /= (real)ns;
	myi1D = 0;
	while (Rmed1D[myi1D] < Rmed[i]) myi1D++;
	myi1D--;
	dens0 = (Rmed1D[myi1D+1]-Rmed[i])*Sigma1D[myi1D] + (Rmed[i]-Rmed1D[myi1D])*Sigma1D[myi1D+1];
	dens0 /= (Rmed1D[myi1D+1]-Rmed1D[myi1D]);
	vrad0 = (Rmed1D[myi1D+1]-Rmed[i])*Vrad1D[myi1D] + (Rmed[i]-Rmed1D[myi1D])*Vrad1D[myi1D+1];
	vrad0 /= (Rmed1D[myi1D+1]-Rmed1D[myi1D]);
	if (DustFluid && (OpenInnerDust == NO)) { /* could be improved... */  
	  /* May 2016: damping towards axisymmetric inst. profiles
	     doesn't to do good for the dust's surface density,
	     changing to initial profiles... */
	  /*
	  dvrad0   = 0.0;
	  dvtheta0 = 0.0;
	  ddens0   = 0.0;
	  for (j = 0; j < ns; j++) {
	    l = i*ns + j;
	    dvrad0   += dvrad[l];
	    dvtheta0 += dvtheta[l];
	    ddens0   += ddens[l];
	  }
	  dvrad0   /= (real)ns;
	  dvtheta0 /= (real)ns;
	  ddens0   /= (real)ns;
	  */
	  dvtheta0 = DVthetaMed[i]-Rmed[i]*OmegaFrame;
	  dvrad0 = DVradMed[i];
	  ddens0 = DSigmaMed[i];
	}
      }
      //
      if (DustFluid && DampToAxi && (OpenInnerDust == NO)) {
	dvrad0   = 0.0;
	dvtheta0 = 0.0;
	ddens0   = 0.0;
	for (j = 0; j < ns; j++) {
	  l = i*ns + j;
	  dvrad0   += dvrad[l];
	  dvtheta0 += dvtheta[l];
	  ddens0   += ddens[l];
	}
	dvrad0   /= (real)ns;
	dvtheta0 /= (real)ns;
	ddens0   /= (real)ns;
      }
      /* Do not modify lines below */
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vrad[l]   = (vrad[l]+lambda*vrad0)/(1.0+lambda);
	vtheta[l] = (vtheta[l]+lambda*vtheta0)/(1.0+lambda);
	dens[l]   = (dens[l]+lambda*dens0)/(1.0+lambda);
	if (EnergyEquation)
	  energ[l]  = (energ[l]+lambda*energ0)/(1.0+lambda);
	if (DustFluid && (OpenInnerDust == NO)) {
	  dvrad[l]   = (dvrad[l]+lambda*dvrad0)/(1.0+lambda);
	  dvtheta[l] = (dvtheta[l]+lambda*dvtheta0)/(1.0+lambda);
	  ddens[l]   = (ddens[l]+lambda*ddens0)/(1.0+lambda);
	}
      }
    }
  }
}

void EvanescentBoundaryDust (DVrad, DVtheta, DRho, step)
     PolarGrid *DVrad, *DVtheta, *DRho;
     real step;
{
  int i, j, l, nr, ns;
  int myi1D;
  real *dvrad, *dvtheta, *ddens;
  real dvrad0, dvtheta0, ddens0;
  real damping, Tin, Tout, lambda;
  dvrad = DVrad->Field;
  dvtheta = DVtheta->Field;
  ddens = DRho->Field;
  nr = DRho->Nrad;
  ns = DRho->Nsec;
  
  lambda = 0.0;
  for (i = 0; i < nr; i++) {
    if ( (Rmed[i] < WKZRMIN) || (Rmed[i] > WKZRMAX) ) {
      /* Damping operates only inside the wave killing zones */
      if (Rmed[i] < WKZRMIN) {
	damping = (Rmed[i]-WKZRMIN)/(GlobalRmed[0]-WKZRMIN);
	Tin = 0.3*pow(Rmed[i],1.5);  // 0.3 Omega_K^-1
	lambda = damping*damping*step/Tin;
      }
      if (Rmed[i] > WKZRMAX) {
	damping = (Rmed[i]-WKZRMAX)/(GlobalRmed[GLOBALNRAD-1]-WKZRMAX);
	Tout = 0.3*pow(Rmed[i],1.5);  // 0.3 Omega_K^-1
	lambda = damping*damping*step/Tout;
      }
      /* Damping wrt initial profiles (requires DampToIni set to yes
	 in .par file) */
      dvrad0   = 0.0;
      dvtheta0 = 0.0;
      ddens0   = 0.0;
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	dvrad0   += dvrad[l];
	dvtheta0 += dvtheta[l];
	ddens0   += ddens[l];
      }
      dvrad0   /= (real)ns;
      dvtheta0 /= (real)ns;
      ddens0   /= (real)ns;
      /* Do not modify lines below */
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	dvrad[l]   = (dvrad[l]+lambda*dvrad0)/(1.0+lambda);
	dvtheta[l] = (dvtheta[l]+lambda*dvtheta0)/(1.0+lambda);
	ddens[l]   = (ddens[l]+lambda*ddens0)/(1.0+lambda);
      }
    }
  }
}

void ApplyOuterSourceMass (Rho, Vrad)
     PolarGrid *Rho, *Vrad;
{
  int i, j, l, nr, ns;
  real *rho, average_rho = 0.0, *vr, penul_vr, average_vr = 0.0;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  if (CPU_Rank == CPU_Highest) {
    i = nr-1;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      average_rho += rho[l];
      average_vr  += vr[l];
    }
    average_rho /= (real)ns;
    average_vr  /= (real)ns;
    average_rho = SigmaMed[nr-1]-average_rho;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rho[l] += average_rho;
    }
    //penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[nr-1]/1.0),SIGMASLOPE-1.0);  
    penul_vr = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[nr-1]/Rinf[nr-1];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      //vr[l] = penul_vr;
      vr[l] += (penul_vr - average_vr);
    }
  }
  if (CPU_Rank == 0) {
    i = 0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      average_rho += rho[l];
      average_vr  += vr[l];
    }
    average_rho /= (real)ns;
    average_vr  /= (real)ns;
    average_rho = SigmaMed[0]-average_rho;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rho[l] += average_rho;
    }
    //penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[0]/1.0),SIGMASLOPE-1.0);  
    penul_vr = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[0]/Rinf[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      //vr[l] = penul_vr;
      vr[l] += (penul_vr - average_vr);
    }
  }
}

void ApplySubKeplerianBoundary (Vtheta, DVtheta)
     PolarGrid *Vtheta, *DVtheta;
     /* New (Jan 2015): we now make use of VthetaMed arrays, defined
	in PframeForce.c */
{
  int i, j, l, nr, ns;
  real *vt, *dvt;
  extern boolean DustFluid;
  vt = Vtheta->Field;
  dvt = DVtheta->Field;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  /* ----- */
  /* i = 0 */
  /* ----- */
  if ( CPU_Rank == 0 ) {
    for (j = 0; j < ns; j++) {
      vt[j] = VthetaMed[0]-Rmed[0]*OmegaFrame;
      if (DustFluid) {
	      dvt[j] = DVthetaMed[0]-Rmed[0]*OmegaFrame;
      }
    }
  }
  /* ---------- */
  /* i = nr - 1 */
  /* ---------- */
  if ( CPU_Rank == CPU_Highest ) {
    i = nr - 1;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      vt[l] = VthetaMed[i]-Rmed[i]*OmegaFrame;
      if (DustFluid) {
	      dvt[l] = DVthetaMed[i]-Rmed[i]*OmegaFrame;
      }
    }
  }
}

void ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, step, sys)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     PolarGrid *DVrad, *DVtheta, *DRho;
     real step;
     PlanetarySystem *sys;
{
  extern boolean DustFluid;
  if (OpenInner == YES) {
    OpenBoundary (Vrad, Vtheta, Rho, Energy);
    if (DustFluid && (OpenInnerDust == YES)) {
      OpenBoundaryd (DVrad, DVtheta, DRho);
      //OpenBoundaryd (DVrad, DVtheta, Rho, DRho);
    }
  }
  // May 2019: outer source mass function shifted prior to call to 
  // non-reflecting BC, otherwise planet wakes bounce off the grid edges
  //if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
  //
  if (NonReflecting == YES) {
    if (EnergyEquation)
      ComputeSoundSpeed (Rho, Energy, sys);
    NonReflectingBoundary (Vrad, Rho, Energy, Vtheta, sys);
    if (DustFluid) {
      if (OpenInnerDust == YES)
	      OpenBoundaryd (DVrad, DVtheta, DRho);
      else
	      EvanescentBoundaryDust (DVrad, DVtheta, DRho, step);
    }
  }
  if (Evanescent == YES) {
    EvanescentBoundary (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, step);
    if (DustFluid && (OpenInnerDust == YES))
      OpenBoundaryd (DVrad, DVtheta, DRho);
  }
  /* New 'mixed' boundary condition, where an open BC is applied at
     the grid's inner edge, and an evanescent BC at the outer edge (WKRMAX
     needs to be specified) */
  if (MixedBC == YES) {
    OpenBoundary (Vrad, Vtheta, Rho, Energy);
    EvanescentBoundary (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, step);
  }
  if (AccBoundary == YES) {
    AccretingBoundary(Vrad, Vtheta, Rho, Energy);
    EvanescentBoundary (Vrad, Vtheta, Rho, Energy, DVrad, DVtheta, DRho, step);
  }
  if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
}

void CorrectVtheta (vtheta, domega)
     PolarGrid *vtheta;
     real domega;
{
  int i, j, l, nr, ns;
  real *vt;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  vt = vtheta->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vt[l] -= domega*Rmed[i];
    }
  }
}

void DampDensity(Vrad, Vtheta, Rho, Energy, step, sys)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
     PlanetarySystem *sys;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real Tadd, axidens;
  real lambda, x0, y0, rp0, x1, y1, rp1;
  real cutrmin, cutrmax, cutrmed, damping, normdampdist;
  real vrad0, vtheta0, energ0, dens0;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* OLD Damping used for adding mass in turbulent runs */
  /*
    for (i = 0; i < nr; i++) {
    // Damping timescale is 20 local orbital periods
    Tadd = 20.0*2.0*M_PI*pow(Rmed[i],1.5);
    axidens = 0.0;
    for (j = 0; j < ns; j++) {
    l = i*ns + j;
    axidens += dens[l];
    }
    // axisymmetric density at ring i 
    axidens /= (real)ns;
    for (j = 0; j < ns; j++) {
    l = i*ns + j;
    dens[l] = dens[l] - (axidens-SigmaMed[i])*dt/Tadd;
    }
    }
  */
  lambda = 0.0;
  x0 = sys->x[0];
  y0 = sys->y[0];
  rp0 = sqrt(x0*x0 + y0*y0);
  x1 = sys->x[1];
  y1 = sys->y[1];
  rp1 = sqrt(x1*x1 + y1*y1);
  cutrmin = rp0*(1.0+DENSDAMPRAD*ASPECTRATIO);
  cutrmax = rp1*(1.0-DENSDAMPRAD*ASPECTRATIO);
  cutrmed = 0.5*(cutrmin+cutrmax);
  for (i = 0; i < nr; i++) {
    if ( (Rmed[i] > cutrmin) && (Rmed[i] < cutrmax) ) {
      if ( (Rmed[i] > cutrmin) && (Rmed[i] < cutrmed) ) 
	normdampdist = (Rmed[i]-cutrmin)/(cutrmed-cutrmin);
      if ( (Rmed[i] >= cutrmed) && (Rmed[i] < cutrmax) )
	normdampdist = (Rmed[i]-cutrmax)/(cutrmed-cutrmax);
      damping = pow(sin((normdampdist)*0.5*M_PI),2.);
      lambda = damping*step;  // to be checked / customized...
      /* damping towards instantaneous axisymmetric fields... */
      vrad0   = 0.0;
      vtheta0 = 0.0;
      dens0   = 0.0;
      energ0  = 0.0;
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vrad0   += vrad[l];
	vtheta0 += vtheta[l];
	dens0   += dens[l];
	energ0  += energ[l];
      }
      vrad0   /= (real)ns;
      vtheta0 /= (real)ns;
      dens0   /= (real)ns;
      energ0  /= (real)ns;
      /* Do not modify lines below */
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vrad[l]   = (vrad[l]+lambda*vrad0)/(1.0+lambda);
	vtheta[l] = (vtheta[l]+lambda*vtheta0)/(1.0+lambda);
	dens[l]   = (dens[l]+lambda*dens0)/(1.0+lambda);
	if (EnergyEquation)
	  energ[l]  = (energ[l]+lambda*energ0)/(1.0+lambda);
      }
    }
  }
}

void Evaporation(Rho, dt)
     PolarGrid *Rho;
     real dt;
{
  int i, j, l, nr, ns;
  real *dens;
  real axidens, floor;
  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  if (ievaporation == 0) {
    /* Calculate axisymmetric density field at beginning of simulation */
    GLOBAL_Axidens_Evap = (real *) malloc(sizeof(real) * GLOBALNRAD);
    if ( GLOBAL_Axidens_Evap==NULL ) {
      fprintf (stderr, "Not enough memory for allocation of GLOBAL_Axidens_Evap in SideEuler.c \n");
      prs_exit (1);
    }
    mpi_make1Dprofile (dens, GLOBAL_Axidens_Evap);
    ievaporation = 1;
  }

  /* Target density = minimum density allowed*/
  for (i = 0; i < nr; i++) {
    floor = 1e-3*GLOBAL_Axidens_Evap[i+IMIN];
    axidens = 0.0;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      axidens += dens[l];
    }
    axidens /= (real)ns;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      dens[l] = dens[l] - (axidens-floor)*dt/TEVAP/2.0/M_PI;
    }
  }
}

/* Function written by S. Ataiee that calculates the scaling factor
 for photevaporation using formula by Owen2012. */
void ApplyPhotoEvaporation (Vrad, Rho, dt)
     PolarGrid *Vrad, *Rho;
     real dt;
{
  int i, j, l, nr, ns;
  int ihole=0;
  real *dens, *vrad, axidens[GLOBALNRAD];
  real x, sigmareduc;
  real scale, Rhole;
  real mdot, summ;
  real world_summ;
  real checkmass;
  real sigcrit;
  real world_mass, world_rhole;
  real l0, lx, l10, coeff;
  real a1,b1,c1,d1,e1,f1,g1;
  real sigd;
  real a2,b2,c2,d2,e2,f2;
  checkmass = 0.0;
  a2 = -0.438226;
  b2 = -0.10658387;
  c2 = 0.5699464;
  d2 = 0.010732277;
  e2 = -0.131809597;
  f2 = -1.32285709;
  l10 = log(10);
  a1 = 0.15138;
  b1 = -1.2182;
  c1 = 3.4046;
  d1 = -3.5717;
  e1 = -0.32762;
  f1 = 3.6064;
  g1 = -2.4918;
  dens = Rho->Field;
  vrad = Vrad->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  /* Calculate axisymmetric density field */
  mpi_make1Dprofile (dens, axidens);

  /* In order to avoid getting stuck before an artificial bump close to 
     the inner boundary, we make the criteria a little fluffy */
  Rhole = GlobalRmed[0];
  if ((axidens[0] <= floordens) || (axidens[1] <= floordens)) 
    sigcrit = 1.1 * floordens;
  else
    sigcrit = floordens;
  while ((axidens[ihole] <= sigcrit) && (ihole < GLOBALNRAD)){
    Rhole = GlobalRmed[ihole];
    ihole++;
  }
  /* ------------------------------ */
  /* Case where hole is in the disc */
  /* ------------------------------ */
  if (Rhole > GlobalRmed[1]) {
    for (i = 0; i <= nr; i++){
      if (Rsup[i] <= Rhole){
	for (j = 0; j <= ns; j++){
	  l = i*ns+j;
	  dens[l] = floordens;
	}
      }
    }
    summ=0;
    for (i = Zero_or_active; i < Max_or_active; i++){
      x = 0.95 * ((Rmed[i]-Rhole)*FACTORUNITLENGTH) * pow(FACTORUNITMASS,-1.0);
      if (x >= 0.0) {
	coeff = a2*b2*exp(b2*x) + c2*d2*exp(d2*x)+ e2*f2*exp(f2*x);
	coeff /= (Rmed[i]*FACTORUNITLENGTH);
	sigd = coeff * exp(-pow((x/57.0),10.0));
      } 
      else
	sigd = 0.0;
      summ += ns*Surf[i]*sigd;  // this is the local Mdot
    }
    MPI_Allreduce (&summ, &world_summ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    summ = world_summ;          // this is the global Mdot summed over all CPUS
    summ *= pow(unit_length*1e2,2);
    mdot = 4.8e-9 * pow(FACTORUNITMASS,-0.148) * pow((LX/1e30),1.14);
    mdot *= 1.9891e33/31556926;  // expected Mdot_wind in g/s
    scale = mdot/summ; 
    scale /= (unit_mass*1e3);    //scale is calculated in cgs
    scale *= pow((unit_length*1e2),2);
    scale *= unit_time; //converting sigmadot to code units       
    for (i = 0; i < nr; i++) {
      x = 0.95 * ((Rmed[i]-Rhole)*FACTORUNITLENGTH) * pow(FACTORUNITMASS,-1.0);
      if (x >= 0.0) {
	coeff = a2*b2*exp(b2*x) + c2*d2*exp(d2*x)+ e2*f2*exp(f2*x);
	coeff /= (Rmed[i]*FACTORUNITLENGTH);
	sigd = coeff * exp(-pow((x/57.0),10.0));
      } 
      else
	sigd = 0.0;
      sigmareduc = scale * sigd ;
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	dens[l] = dens[l] - sigmareduc*dt;
	if ((dens[l] <= floordens)){
	  dens[l] = floordens;
	} else {
	  checkmass += sigmareduc * Surf[i] * dt;
	}
      }
    }
  } else {
    /* ----------------------------------- */
    /* Case where hole is outside the grid */
    /* ----------------------------------- */
    summ = 0;
    for (i = Zero_or_active; i < Max_or_active; i++){
      x = 0.85 * (Rmed[i]*FACTORUNITLENGTH) * pow(FACTORUNITMASS,-1.0);
      l0 = log10(x);
      lx = log(x);
      if (x > 0.7){
	sigd = pow(10.0,(a1*pow(l0,6.0) + b1*pow(l0,5.0) + c1*pow(l0,4.0)+ d1*pow(l0,3.0) \
		       + e1*pow(l0,2.0) + f1*l0 + g1));
	coeff = 6.0*a1*pow(lx,5.0)/pow(x,2.0)/pow(l10,7.0); 
	coeff+= 5.0*b1*pow(lx,4.0)/pow(x,2.0)/pow(l10,6.0);
	coeff+= 4.0*c1*pow(lx,3.0)/pow(x,2.0)/pow(l10,5.0);
	coeff+= 3.0*d1*pow(lx,2.0)/pow(x,2.0)/pow(l10,4.0);
	coeff+= 2.0*e1*lx/pow(x,2.0)/pow(l10,3.0);
	coeff+= 1.0*f1/pow(x,2.0)/pow(l10,2.0);
	sigd *= (coeff*exp(-pow((x/100.),10.0)));
      } 
      else
	sigd = 0.0;
      summ += ns*Surf[i]*sigd;
    }
    MPI_Allreduce (&summ, &world_summ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    summ = world_summ;
    summ *= pow(unit_length*1e2,2);
    mdot = 6.25e-9 * pow(FACTORUNITMASS,-0.068) * pow((LX/1e30),1.14);
    mdot *= 1.9891e33/31556926;
    scale = mdot/summ;
    scale /= (unit_mass*1e3); //scale is calculated in cgs
    scale *= pow((unit_length*1e2),2);
    scale *= unit_time; //converting sigmadot to code units   
    for (i = 0; i < nr; i++) {
      x = 0.85 * (Rmed[i]*FACTORUNITLENGTH) * pow(FACTORUNITMASS,-1.0);
      l0 = log10(x);
      lx = log(x);
      if ( x > 0.7 ){
	sigd = pow(10.0,(a1*pow(l0,6.0) + b1*pow(l0,5.0) + c1*pow(l0,4.0)+ d1*pow(l0,3.0) \
		       + e1*pow(l0,2.0) + f1*l0 + g1));
	coeff = 6.0*a1*pow(lx,5.0)/pow(x,2.0)/pow(l10,7.0); 
	coeff+= 5.0*b1*pow(lx,4.0)/pow(x,2.0)/pow(l10,6.0);
	coeff+= 4.0*c1*pow(lx,3.0)/pow(x,2.0)/pow(l10,5.0);
	coeff+= 3.0*d1*pow(lx,2.0)/pow(x,2.0)/pow(l10,4.0);
	coeff+= 2.0*e1*lx/pow(x,2.0)/pow(l10,3.0);
	coeff+= 1.0*f1/pow(x,2.0)/pow(l10,2.0);
	sigd *= (coeff*exp(-pow((x/100.),10.0)));
      } 
      else
	sigd = 0;
      sigmareduc = scale * sigd ;
      //printf ("CPU %d: at r=%lg, dsigma/sigma=%lg\n",CPU_Rank,Rmed[i],sigmareduc*dt/axidens[i+IMIN]);
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	dens[l] = dens[l] - sigmareduc*dt;
	if (dens[l] <= floordens){
	  dens[l] = floordens;
	} else {
	  checkmass += sigmareduc * Surf[i] * dt;
	}
      }
    }
  }
  MPI_Allreduce (&checkmass, &world_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  checkmass = world_mass;
}


void InitializeOneDViscousEvolution ()
{
  FILE *fich;
  char file[200];
  int i;

  extern boolean Restart;
  extern int     NbRestart;

  Radii1D = (double *) malloc(sizeof(double)*(NRAD1D+1));
  Rmed1D = (double *) malloc(sizeof(double)*NRAD1D);
  Rsup1D = (double *) malloc(sizeof(double)*NRAD1D);
  Rinf1D = (double *) malloc(sizeof(double)*NRAD1D);
  Sigma1D = (double *) malloc(sizeof(double)*NRAD1D);
  Vrad1D = (double *) malloc(sizeof(double)*NRAD1D);
  cs1D = (double *) malloc(sizeof(double)*NRAD1D);
  viscosity1D = (double *) malloc(sizeof(double)*NRAD1D);
  term1 = (double *) malloc(sizeof(double)*NRAD1D);
  term2 = (double *) malloc(sizeof(double)*NRAD1D);
  term3 = (double *) malloc(sizeof(double)*NRAD1D);
  
  for (i=0; i<=NRAD1D; i++) {
    // by default, a logarithmic radial spacing is assumed for the 1D grid
    Radii1D[i] = RMIN1D*exp((real)i/(real)NRAD1D*log(RMAX1D/RMIN1D));
  }
  
  for (i=0; i<NRAD1D; i++) {
    Rinf1D[i] = Radii1D[i];
    Rsup1D[i] = Radii1D[i+1];
    Rmed1D[i] = 2.0/3.0*(Radii1D[i+1]*Radii1D[i+1]*Radii1D[i+1]-Radii1D[i]*Radii1D[i]*Radii1D[i]);
    Rmed1D[i] = Rmed1D[i] / (Radii1D[i+1]*Radii1D[i+1]-Radii1D[i]*Radii1D[i]);
    Sigma1D[i] = Sigma(Rmed1D[i]);
    cs1D[i] = AspectRatio(Rmed1D[i])/sqrt(Rmed1D[i])*pow(Rmed1D[i], FLARINGINDEX);
    if (VISCOSITY != 0.0)
      viscosity1D[i] = VISCOSITY;
    if (ViscosityAlpha) {
      viscosity1D[i] = ALPHAVISCOSITY*cs1D[i]*cs1D[i]*pow(Rmed1D[i],1.5);
      Vrad1D[i] = -3.0*viscosity1D[i]/Rmed1D[i]*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
    }
    else
      Vrad1D[i] = -3.0*viscosity1D[i]/Rmed1D[i]*(-SIGMASLOPE+0.5);
    term1[i] = term2[i] = term3[i] = 0.0;
  }

  /* Restart case */
  if (Restart) {
    sprintf (file, "%s/gasdens.ascii_rad.%d.dat", OUTPUTDIR, NbRestart);
    fich = fopen (file, "r");
    for ( i = 0 ; i < NRAD1D; i++ ) {
      fscanf (fich, "%lg\n%lg", &Rmed1D[i], &Sigma1D[i]);
    }
    fclose(fich);
    sprintf (file, "%s/gasvrad.ascii_rad.%d.dat", OUTPUTDIR, NbRestart);
    fich = fopen (file, "w");
    for ( i = 0 ; i < NRAD1D; i++ ) {
      fscanf (fich, "%lg\n%lg", &Rmed1D[i], &Vrad1D[i]);
    }
    fclose(fich);
  }
  
  if (CPU_Master) {
    sprintf (file, "%s/dims1D.dat", OUTPUTDIR);
    fich = fopen (file, "w");
    fprintf(fich,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",NRAD1D,NRAD1D,NRAD1D,NRAD1D,NRAD1D,NRAD1D,NRAD1D,NRAD1D);
    fclose(fich);
    sprintf (file, "%s/used_rad1D.dat", OUTPUTDIR);
    fich = fopen (file, "w");
    for (i=0; i<=NRAD1D; i++) {
      fprintf(fich,"%f\n",Radii1D[i]);
    }
    fclose(fich);
  }
}

void SolveOneDViscousEvolution (dt)
     real dt;
{
  int i;
  real drm, drp, drnu, drsigma, d2rnu, d2rsigma;
  extern boolean BC1D_SS_ZeroVel, BC1D_SS_NonZeroVel, BC1D_ZeroDens;

  /* Update of gas surface density and radial velocity. Note that viscosity1D isn't
     evolved in time, it is set by the isothermal sound speed via the
     parameters in the .par file */
  for (i=1; i<=NRAD1D-2; i++) {
    drm = Rmed1D[i]-Rmed1D[i-1];
    drp = Rmed1D[i+1]-Rmed1D[i];
    drnu     = (drm*drm*viscosity1D[i+1]-drp*drp*viscosity1D[i-1] + (drp*drp-drm*drm)*viscosity1D[i]) / (drm*drm*drp + drp*drp*drm);
    drsigma  = (drm*drm*Sigma1D[i+1]-drp*drp*Sigma1D[i-1] + (drp*drp-drm*drm)*Sigma1D[i]) / (drm*drm*drp + drp*drp*drm);
    d2rnu    = (drm*viscosity1D[i+1] + drp*viscosity1D[i-1] - (drm+drp)*viscosity1D[i]) / (0.5*drm*drp*drp + 0.5*drp*drm*drm); 
    d2rsigma = (drm*Sigma1D[i+1] + drp*Sigma1D[i-1] - (drm+drp)*Sigma1D[i]) / (0.5*drm*drp*drp + 0.5*drp*drm*drm); 
    term1[i] = 3.0*Sigma1D[i]*(1.5*drnu/Rmed1D[i] + d2rnu);
    term2[i] = 3.0*drsigma*(1.5*viscosity1D[i]/Rmed1D[i] + 2.0*drnu);
    term3[i] = 3.0*viscosity1D[i]*d2rsigma;
    Vrad1D[i] = -1.5*viscosity1D[i]/Rmed1D[i] - 3.0*drnu - 1.5*(viscosity1D[i]/Sigma1D[i])*drsigma;
  }
  for (i=1; i<=NRAD1D-2; i++) {
    Sigma1D[i] += dt*(term1[i] + term2[i] + term3[i]);
  }
 
  /* Boundary conditions */
  if (BC1D_SS_ZeroVel) {
    // Apply zero radial velocity steady state condition
    Sigma1D[0] = Sigma1D[1] * (viscosity1D[1]/viscosity1D[0]) * sqrt(Rmed1D[1]/Rmed1D[0]);
    Sigma1D[NRAD1D-1] = Sigma1D[NRAD1D-2] * (viscosity1D[NRAD1D-2]/viscosity1D[NRAD1D-1]) * sqrt(Rmed1D[NRAD1D-2]/Rmed1D[NRAD1D-1]);
  }
  if (BC1D_SS_NonZeroVel) {
    // Apply non-zero radial velocity steady state condition vr = -3/2 nu/R
    Sigma1D[0] = Sigma1D[1] * (viscosity1D[1]/viscosity1D[0]);
    Sigma1D[NRAD1D-1] = Sigma1D[NRAD1D-2] * (viscosity1D[NRAD1D-2]/viscosity1D[NRAD1D-1]);
  }
  if (BC1D_ZeroDens) {
    // Apply zero-torque 'low-density' boundary condition
    Sigma1D[0] = 1e-9;
    Sigma1D[NRAD1D-1] = 1e-9;
  }
  Vrad1D[0] = Vrad1D[1];
  Vrad1D[NRAD1D-1] = Vrad1D[NRAD1D-2];
}

void Write1DViscProfiles (number)
     int number;
{
  FILE *fich;
  char file[200];
  int i;
  if (CPU_Master) {
    sprintf (file, "%s/gasdens.ascii_rad.%d.dat", OUTPUTDIR, number);
    fich = fopen (file, "w");
    for ( i = 0 ; i < NRAD1D; i++ ) {
      fprintf(fich,"%lg\t%lg\n",Rmed1D[i], Sigma1D[i]);
    }
    fclose(fich);
    sprintf (file, "%s/gasvrad.ascii_rad.%d.dat", OUTPUTDIR, number);
    fich = fopen (file, "w");
    for ( i = 0 ; i < NRAD1D; i++ ) {
      fprintf(fich,"%lg\t%lg\n",Rmed1D[i], Vrad1D[i]);
    }
    fclose(fich);
  }
}
