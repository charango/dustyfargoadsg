/** \file Viscosity.c

Calculation of the viscous %force.  The function FViscosity() returns
the (kinematic) viscosity as a function of the radius (it handles all
case: alpha or uniform viscosity, and inner cavity with a different
viscosity). The calculation of the stress tensor elements in 2D
cylindrical coordinates is done in ComputeViscousTerms(). The update
of the velocity is done in UpdateVelocitiesWithViscosity(). This file
also contains the function AspectRatio(), which gives the aspect ratio
as a function of the radius, in the case of a temperature jump in the
disk (much in the manner as cavities arising from a viscosity jump are
handled, hence the location of this function). Note that AspectRatio()
does not feature the FLARINGINDEX, which is taken into account by the
calling function.

*/

#include "mp.h"

static PolarGrid *DRR, *DRP, *DPP;
static PolarGrid *DDRR, *DDRP, *DDPP;

real FViscosity (rad)
     real rad;
{
  real viscosity, rmin, rmax, scale, nudot, nuint, nuext, fdrop;
  extern boolean TailOffSareh;
  int i = 0;
  viscosity = VISCOSITY;
  if (ViscosityAlpha) {
    while (GlobalRmed[i] < rad) i++;
    /* GLOBAL_SoundSpeed contains global, axisymmetric soundspeed
       array */
    viscosity = ALPHAVISCOSITY*GLOBAL_SoundSpeed[i]*	\
      GLOBAL_SoundSpeed[i]*pow(rad, 1.5);
  }
  rmin = CAVITYRADIUS-CAVITYWIDTH*ASPECTRATIO*pow(CAVITYRADIUS,1.0+FLARINGINDEX);
  rmax = CAVITYRADIUS+CAVITYWIDTH*ASPECTRATIO*pow(CAVITYRADIUS,1.0+FLARINGINDEX);
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) viscosity *= CAVITYRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    viscosity *= exp((rmax-rad)/(rmax-rmin)*log(CAVITYRATIO));
    /*
    # viscosity profile used in Debras+ 21
    nuint = CAVITYRATIO*VISCOSITY;
    nuext = VISCOSITY;
    viscosity = nuint + (nuext-nuint)*log(rad/rmin)/log(rmax/rmin);
    */
  }
  if (TailOffSareh) {
    fdrop = 1.0/(1.0 + exp(-(rad-1.5)/0.1));
    viscosity /= (fdrop + 0.01*(1.0-fdrop));
  }

  /* Linear time variation of viscosity over time duration
     RELEASEDATEVISCOSITY between VISCOSITY and RELEASEVISCOSITY */
  if ( (!ViscosityAlpha) && (RELEASEDATEVISCOSITY > 1e-3) ) {
    if (PhysicalTime < RELEASEDATEVISCOSITY) {
      nudot = (RELEASEVISCOSITY-VISCOSITY) / (RELEASEDATEVISCOSITY-PhysicalTimeInitial);
      viscosity = VISCOSITY + nudot*PhysicalTime;
    } else {
      viscosity = RELEASEVISCOSITY;
    }
  }
  return viscosity;
}

real DFViscosity (rad)
     real rad;
{
  real viscosity, rmin, rmax, scale, nudot;
  int i = 0;
  viscosity = DVISCOSITY;
  if (DViscosityAlpha) {
    while (GlobalRmed[i] < rad) i++;
    /* GLOBAL_DustSoundSpeed contains global, axisymmetric soundspeed
       array for the dust */
    viscosity = DALPHAVISCOSITY*GLOBAL_DustSoundSpeed[i]*	\
      GLOBAL_DustSoundSpeed[i]*pow(rad, 1.5);
  }
  rmin = CAVITYRADIUS-CAVITYWIDTH*DASPECTRATIO;
  rmax = CAVITYRADIUS+CAVITYWIDTH*DASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) viscosity *= CAVITYRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    viscosity *= exp((rmax-rad)/(rmax-rmin)*log(CAVITYRATIO));
  }
  /* Linear time variation of viscosity over time duration
     RELEASEDATEVISCOSITY between DVISCOSITY and RELEASEVISCOSITY */
  if (RELEASEDATEVISCOSITY > 1e-3) {
    if (PhysicalTime < RELEASEDATEVISCOSITY) {
      nudot = (RELEASEVISCOSITY-DVISCOSITY) / (RELEASEDATEVISCOSITY-PhysicalTimeInitial);
      viscosity = DVISCOSITY + nudot*PhysicalTime;
    } else {
      viscosity = RELEASEVISCOSITY;
    }
  }
  return viscosity;
}


real AspectRatio (rad)
     real rad;
{
  real aspectratio, rmin, rmax, scale;
  aspectratio = ASPECTRATIO;
  rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*ASPECTRATIO;
  rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*ASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) aspectratio *= TRANSITIONRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    aspectratio *= exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
  }
  return aspectratio;
}

real DAspectRatio (rad)
     real rad;
{
  real aspectratio, rmin, rmax, scale;
  aspectratio = DASPECTRATIO;
  rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*DASPECTRATIO;
  rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*DASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) aspectratio *= TRANSITIONRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    aspectratio *= exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
  }
  return aspectratio;
}

void InitViscosity ()
{
  DivergenceVelocity = CreatePolarGrid(NRAD, NSEC, "DivV");
  DRR                = CreatePolarGrid(NRAD, NSEC, "Drr");
  DRP                = CreatePolarGrid(NRAD, NSEC, "Drp");
  DPP                = CreatePolarGrid(NRAD, NSEC, "Dpp");
  TAURR              = CreatePolarGrid(NRAD, NSEC, "TAUrr");
  TAURP              = CreatePolarGrid(NRAD, NSEC, "TAUrp");
  TAUPP              = CreatePolarGrid(NRAD, NSEC, "TAUpp");
  DDivergenceVelocity = CreatePolarGrid(NRAD, NSEC, "DDivV");
  DDRR                = CreatePolarGrid(NRAD, NSEC, "DDrr");
  DDRP                = CreatePolarGrid(NRAD, NSEC, "DDrp");
  DDPP                = CreatePolarGrid(NRAD, NSEC, "DDpp");
  DTAURR              = CreatePolarGrid(NRAD, NSEC, "DTAUrr");
  DTAURP              = CreatePolarGrid(NRAD, NSEC, "DTAUrp");
  DTAUPP              = CreatePolarGrid(NRAD, NSEC, "DTAUpp");
}

void ComputeViscousTerms (RadialVelocity, AzimuthalVelocity, Rho, DRadialVelocity, DAzimuthalVelocity, DRho)
     PolarGrid *RadialVelocity, *AzimuthalVelocity, *Rho;
     PolarGrid *DRadialVelocity, *DAzimuthalVelocity, *DRho;
{
  int i, j, l, nr, ns;
  int lip, ljp, ljm, lim, ljmim;
  real *rho, *vr, *vt, *cs;
  real *Drr, *Drp, *Dpp, *divergence;
  real *Trr, *Trp, *Tpp;
  real *Drho, *Dvr, *Dvt, *Dcs;
  real *DDrr, *DDrp, *DDpp, *Ddivergence;
  real *DTrr, *DTrp, *DTpp;
  real dphi, onethird, invdphi;
  real viscosity, Dviscosity;
  extern boolean DustFluid;
  ComputeDivergenceVelocity (RadialVelocity, AzimuthalVelocity, DRadialVelocity, DAzimuthalVelocity);
  nr  = Rho->Nrad;
  ns  = Rho->Nsec;
  rho = Rho->Field;
  vr  = RadialVelocity->Field;
  vt  = AzimuthalVelocity->Field;
  divergence = DivergenceVelocity->Field;
  Drr = DRR->Field;
  Drp = DRP->Field;
  Dpp = DPP->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  cs = SoundSpeed->Field;
  // Dust Fluid
  if (DustFluid) {
    Drho = DRho->Field;
    Dvr  = DRadialVelocity->Field;
    Dvt  = DAzimuthalVelocity->Field;
    Ddivergence = DDivergenceVelocity->Field;
    DDrr = DDRR->Field;
    DDrp = DDRP->Field;
    DDpp = DDPP->Field;
    DTrr = DTAURR->Field;
    DTrp = DTAURP->Field;
    DTpp = DTAUPP->Field;
    Dcs = DSoundSpeed->Field;
  }
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  onethird = 1.0/3.0;
  if (ViscosityAlpha)
    mpi_make1Dprofile (cs, GLOBAL_SoundSpeed);
  if (DustFluid && DViscosityAlpha)
    mpi_make1Dprofile (Dcs, GLOBAL_DustSoundSpeed);
#pragma omp parallel private(l,lip,ljp,j,ljm,lim)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* Drr, Dpp and divV computation */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lip = l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	Drr[l] = (vr[lip]-vr[l])*InvDiffRsup[i];
	Dpp[l] = (vt[ljp]-vt[l])*invdphi*InvRmed[i]+0.5*(vr[lip]+vr[l])*InvRmed[i];
	// Dust Fluid
	if (DustFluid) {
	  DDrr[l] = (Dvr[lip]-Dvr[l])*InvDiffRsup[i];
	  DDpp[l] = (Dvt[ljp]-Dvt[l])*invdphi*InvRmed[i]+0.5*(Dvr[lip]+Dvr[l])*InvRmed[i];
	}
      }
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* Drp computation */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	lim = l-ns;
	Drp[l] = 0.5*(Rinf[i]*(vt[l]*InvRmed[i]-vt[lim]*InvRmed[i-1])*InvDiffRmed[i]+ \
		      (vr[l]-vr[ljm])*invdphi*InvRinf[i]);
	// Dust Fluid
	if (DustFluid) {
	  DDrp[l] = 0.5*(Rinf[i]*(Dvt[l]*InvRmed[i]-Dvt[lim]*InvRmed[i-1])*InvDiffRmed[i]+ \
			 (Dvr[l]-Dvr[ljm])*invdphi*InvRinf[i]);
	}
      }
    }
  }
#pragma omp parallel private(l,ljmim,j,ljm,lim,viscosity)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* TAUrr and TAUpp computation */
      if ((ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0))
	viscosity = FViscosity(Rmed[i]);
      else
	viscosity = 0.0;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	// divergence term is now calculated in ComputeDivergenceVelocity function
	Trr[l] = 2.0*rho[l]*viscosity*(Drr[l]-onethird*divergence[l]);
	Tpp[l] = 2.0*rho[l]*viscosity*(Dpp[l]-onethird*divergence[l]);
      }
      // Dust Fluid
      if (DustFluid) {
	if ( (DALPHAVISCOSITY != 0.0) || (DVISCOSITY != 0.0) )
	  Dviscosity = DFViscosity (Rmed[i]);
	else
	  Dviscosity = 0.0;
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  // divergence term is now calculated in ComputeDivergenceVelocity function
	  DTrr[l] = 2.0*Drho[l]*Dviscosity*(DDrr[l]-onethird*Ddivergence[l]);
	  DTpp[l] = 2.0*Drho[l]*Dviscosity*(DDpp[l]-onethird*Ddivergence[l]);
	}
      }
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* TAUrp computation */
      if ((ALPHAVISCOSITY != 0.0) || (VISCOSITY != 0.0))
	viscosity = FViscosity (Rmed[i]);
      else
	viscosity = 0.0;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljmim=ljm-ns;
	Trp[l] = 2.0*0.25*(rho[l]+rho[lim]+rho[ljm]+rho[ljmim])*viscosity*Drp[l];
      }
      // Dust Fluid
      if (DustFluid) {
	if ( (DALPHAVISCOSITY != 0.0) || (DVISCOSITY != 0.0) )
	  Dviscosity = DFViscosity (Rmed[i]);
	else
	  Dviscosity = 0.0;
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  lim = l-ns;
	  ljm = l-1;
	  if (j == 0) ljm = i*ns+ns-1;
	  ljmim=ljm-ns;
	  DTrp[l] = 2.0*0.25*(Drho[l]+Drho[lim]+Drho[ljm]+Drho[ljmim])*Dviscosity*DDrp[l];
	}
      }
    }
  }
}

void UpdateVelocitiesWithViscosity (RadialVelocity, AzimuthalVelocity, Rho, DRadialVelocity, DAzimuthalVelocity, DRho, DeltaT)
     PolarGrid *RadialVelocity, *AzimuthalVelocity, *Rho;
     PolarGrid *DRadialVelocity, *DAzimuthalVelocity, *DRho;
     real DeltaT;
{
  int i, j, l, nr, ns;
  int lip, ljp, ljm, lim;
  real *rho, *vr, *vt;
  real *Trr, *Trp, *Tpp;
  real *Drho, *Dvr, *Dvt;
  real *DTrr, *DTrp, *DTpp;
  real dphi, invdphi;
  extern boolean DustFluid;
  nr  = Rho->Nrad;
  ns  = Rho->Nsec;
  rho = Rho->Field;
  vr  = RadialVelocity->Field;
  vt  = AzimuthalVelocity->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  if (DustFluid) {
    Drho = DRho->Field;
    Dvr  = DRadialVelocity->Field;
    Dvt  = DAzimuthalVelocity->Field;
    DTrr = DTAURR->Field;
    DTrp = DTAURP->Field;
    DTpp = DTAUPP->Field;
  }
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  /* Now we can update velocities */
  /* with the viscous source term */
  /* of Navier-Stokes equations */
#pragma omp parallel private(l,j,lip,ljp,ljm,lim)
  {
#pragma omp for nowait
    for (i = 1; i < nr-1; i++) {	/* vtheta first */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lip = l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vt[l] += DeltaT*InvRmed[i]*((Rsup[i]*Trp[lip]-Rinf[i]*Trp[l])*InvDiffRsup[i]+ \
				    (Tpp[l]-Tpp[ljm])*invdphi+		\
				    0.5*(Trp[l]+Trp[lip]))/(0.5*(rho[l]+rho[ljm]));
	if (DustFluid) {
	  Dvt[l] += DeltaT*InvRmed[i]*((Rsup[i]*DTrp[lip]-Rinf[i]*DTrp[l])*InvDiffRsup[i]+ \
				       (DTpp[l]-DTpp[ljm])*invdphi+		\
				       0.5*(DTrp[l]+DTrp[lip]))/(0.5*(Drho[l]+Drho[ljm]));
	}
      }
    }
#pragma omp for nowait
    for (i = 1; i < nr; i++) {	/* and now vrad */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vr[l] += DeltaT*InvRinf[i]*((Rmed[i]*Trr[l]-Rmed[i-1]*Trr[lim])*InvDiffRmed[i]+ \
				    (Trp[ljp]-Trp[l])*invdphi-		\
				    0.5*(Tpp[l]+Tpp[lim]))/(0.5*(rho[l]+rho[lim]));
	if (DustFluid) {
	  Dvr[l] += DeltaT*InvRinf[i]*((Rmed[i]*DTrr[l]-Rmed[i-1]*DTrr[lim])*InvDiffRmed[i]+ \
				       (DTrp[ljp]-DTrp[l])*invdphi-		\
				       0.5*(DTpp[l]+DTpp[lim]))/(0.5*(Drho[l]+Drho[lim]));
	}
      }
    }
  }
}
