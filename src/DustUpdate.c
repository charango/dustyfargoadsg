/* Contains leapfrog integrator routines: update of particles
   positions and velocities separately */

#include "mp.h"

void SemiUpdateDustPositions(dsys, psys, gasdens, timestep)
     DustSystem *dsys;
     PlanetarySystem *psys;
     PolarGrid *gasdens;
     real timestep;
{
  int k, i, myk, indexprev=0, indexnext = 0;
  extern boolean DustFeelTurb, ZZIntegrator, RemoveDustFromPlanetsHillRadius;
  real r, visc, dr, psi;
  real x1, x2, u, w;
  real xp, yp, rp, tp, mp, rh, rd, td, distfromplanet, xd, yd;
  real *senddusttoprevCPU, *senddusttonextCPU, *recvdustfromprevCPU, *recvdustfromnextCPU;
  real *dens;
  MPI_Request req1, req2, req3, req4;
  FILE *file;
  char name[256];
  real D_dust, dr_mean, sigma_r_squared, dphi_mean, sigma_phi_squared, dxtheta, invdxtheta, dphi;
  int ns, ip, jp, lp, lpm1, lpjm;
  senddusttoprevCPU = dsys->senddusttoprevCPU;
  senddusttonextCPU = dsys->senddusttonextCPU;
  recvdustfromprevCPU = dsys->recvdustfromprevCPU;
  recvdustfromnextCPU = dsys->recvdustfromnextCPU;
  dens = gasdens->Field;
  ns = gasdens->Nsec;
  for (k=0; k<NBPART; k++) {
    r = dsys->r[k];    //  keep track of R(t) for azimuth's update
    if ( (r >= Rinf[Zero_or_active]) && (r < Rsup[Max_or_active-1]) ) {
      dsys->r[k] += dsys->vr[k]*0.5*timestep;
      /* Zhaohuan Zhu's integrator uses l=R vphi instead of vphi, but both lines are 
	 strictly equivalent */
      if (!ZZIntegrator)
	dsys->th[k] += dsys->vth[k]*r*0.25*timestep*(pow(r,-2.) + pow(dsys->r[k],-2.));
      else
	dsys->th[k] += dsys->l[k]*0.25*timestep*(pow(r,-2.) + pow(dsys->r[k],-2.));
      /* Retain particle's azimuths between AziInf[0] and
	 AziSup[ns-1], the first and last azimuthal interfaces of the
	 grid cells. With a 2PI azimuthal extent, AziInf[0] = 0 -
	 pi/ns, and AziSup[ns-1] = 2pi-pi/ns. */
      if (dsys->th[k] < AziInf[0])
	dsys->th[k] += (AziSup[NSEC-1]-AziInf[0]);  // add 2pi
      if (dsys->th[k] > AziSup[NSEC-1])
	dsys->th[k] -= (AziSup[NSEC-1]-AziInf[0]);  // subtract 2pi
      if ( (dsys->th[k] < AziInf[0]) || (dsys->th[k] > AziSup[NSEC-1]) )
	printf ("Pb in SemiUpdateDustPositions for particle %d: old radius=%lg, new radius=%lg, azimuth=%lg, vr=%lg, vphi=%lg\n",k,r,dsys->r[k],dsys->th[k],dsys->vr[k],dsys->vth[k]);
      /* ----------------------------------------------------------- */
      /* NEW Feb. 2015: add stochastic turbulence to particles, using
	 Charnoz+ prescription. We apply such prescription only if
	 particle's new radius is in current CPU! */
      /* ----------------------------------------------------------- */
      // 04/2020 not sure why I require r to be bound in active CPU...
      if (DustFeelTurb && (dsys->r[k] >= Rinf[Zero_or_active]) && (dsys->r[k] <= Rsup[Max_or_active-1]) ) {
	visc = FViscosity(dsys->r[k]);
	/*
	if ( (dsys->r[k] >= Rinf[Zero_or_active]) && (dsys->r[k] <= Rsup[Max_or_active-1]) )
	  visc = FViscosity(dsys->r[k]);
	else {
	  if (dsys->r[k] < Rinf[Zero_or_active])
	    visc = FViscosity( Rinf[Zero_or_active]);
	  if (dsys->r[k] > Rsup[Max_or_active-1])
	    visc = FViscosity(Rsup[Max_or_active-1]);
	}
	*/
	D_dust  = visc*(1.0+4.0*pow(dsys->stokesnb[k],2.))*pow(1.+pow(dsys->stokesnb[k],2.),-2.0); // dust diffusion coefficient
	ip = Zero_or_active;
	// CB: added second condition Nov. 2016 after segmentation fault has been found when rp > r_outer_edge...
	// ZZ added third condition Nov. 2017 after segmentation fault has been found when rp < Rinf[Zero_or_active]...
	if (dsys->r[k] >= Rinf[ip]) {
	  while ( (dsys->r[k] >= Rinf[ip]) && (ip <= Max_or_active-1) ) ip++;
	}
	ip--;
	if (ip < 0) 
	  printf ("Beware: ip is negative for particle %d, for which ip=%d\n",k,ip);
	tp = dsys->th[k];
	jp = floor(ns*(tp-AziInf[0])/(AziSup[ns-1]-AziInf[0]));
	if (jp == ns)
	  jp = 0;
	lp = ip*ns+jp;
	if (ip != 0)
	  lpm1 = lp-ns;
	else
	  lpm1 = lp;
	lpjm = lp-1;
	if (jp == 0) 
	  lpjm = ip*ns+ns-1;
	// Eq. (21) of Charnoz et al. (2011) applied to r-direction, dD_dust/dr discarded
	dr_mean = D_dust/dens[lp] * (dens[lp]-dens[lpm1])*InvDiffRmed[ip] * timestep;
	sigma_r_squared = 2.0*D_dust*timestep;
	/* dr is sorted with a gaussian distribution of mean
	   dr_mean and standard deviation sigma_r_squared */
	dr = myrandn_rad(dr_mean,sqrt(sigma_r_squared));
	dsys->r[k] = dsys->r[k] + dr;
	// Eq. (21) of Charnoz et al. (2011) applied to phi-direction, dD_dust/dphi discarded
	dxtheta = (PMAX-PMIN)/(real)ns*Rmed[ip];
	invdxtheta = 1.0/dxtheta;
	/* dphi is sorted with a gaussian distribution of mean
	   dphi_mean and standard deviation sigma_r_squared/r */
	dphi_mean = D_dust/dens[lp] * (dens[lp]-dens[lpjm])*invdxtheta * timestep /dsys->r[k];
	dphi = myrandn_phi(dphi_mean,sqrt(sigma_r_squared)/dsys->r[k]);
	dsys->th[k] = dsys->th[k] + dphi;
      }
      /* ----------------------------------------------------------- */
      /* NEW (April 2016): remove particles inside the planet's Hill
	 radius */
      if (RemoveDustFromPlanetsHillRadius) {
	for (i=0; i<psys->nb; i++) {
	  xp = psys->x[i];
	  yp = psys->y[i];
	  rp = sqrt( xp*xp + yp*yp );
	  tp = atan2(yp,xp);
	  distfromplanet = sqrt( dsys->r[k]*dsys->r[k] + rp*rp - 2.0*dsys->r[k]*rp*cos(dsys->th[k]-tp));
	  mp = psys->mass[i];
	  rh = rp*pow(mp/3.0,1./3.);
	  while (distfromplanet < rh) {
	    dsys->th[k] = AziInf[0] + drand48()*(AziSup[NSEC-1]-AziInf[0]);
	    dsys->r[k] = Rinf[Zero_or_active] + drand48()*(Rsup[Max_or_active-1]-Rinf[Zero_or_active]);
	    distfromplanet = sqrt( dsys->r[k]*dsys->r[k] + rp*rp - 2.0*dsys->r[k]*rp*cos(dsys->th[k]-tp));
	  }
	}
      }

      /* Again, retain particle's azimuths between AziInf[0] and
	 AziSup[ns-1], the first and last azimuthal interfaces of the
	 grid cells. With a 2PI azimuthal extent, AziInf[0] = 0 -
	 pi/ns, and AziSup[ns-1] = 2pi-pi/ns. */
      if (dsys->th[k] < AziInf[0])
	dsys->th[k] += (AziSup[NSEC-1]-AziInf[0]);  // add 2pi
      if (dsys->th[k] > AziSup[NSEC-1])
	dsys->th[k] -= (AziSup[NSEC-1]-AziInf[0]);  // subtract 2pi
      
      if ( (dsys->th[k] < AziInf[0]) || (dsys->th[k] > AziSup[NSEC-1]) )
	printf ("Pb in SemiUpdateDustPositions just after DustTurbulence and RemoveDustFromHillRadius for particle %d: old radius=%lg, new radius=%lg, azimuth=%lg, vr=%lg, vphi=%lg, dphi=%lg, dphi_mean=%lg\n",k,r,dsys->r[k],dsys->th[k],dsys->vr[k],dsys->vth[k], dphi, dphi_mean);
      
      
      /* Here we check if the particle's new radius makes it sweep to
	 another CPU. No need to communicate particles size as long as
	 it doesn't change! */
      if ( (dsys->r[k] <  Rinf[Zero_or_active]) && (CPU_Rank != 0) ) {
	senddusttoprevCPU[5*indexprev] = (real)k;
	senddusttoprevCPU[5*indexprev+1] = dsys->r[k];
	senddusttoprevCPU[5*indexprev+2] = dsys->th[k];
	senddusttoprevCPU[5*indexprev+3] = dsys->vr[k];
	if (!ZZIntegrator)
	  senddusttoprevCPU[5*indexprev+4] = dsys->vth[k];
	else
	  senddusttoprevCPU[5*indexprev+4] = dsys->l[k];
	indexprev+=1;
      }
      if ( (dsys->r[k] >= Rsup[Max_or_active-1]) && (CPU_Rank != CPU_Highest) ) {
	senddusttonextCPU[5*indexnext] = (real)k;
	senddusttonextCPU[5*indexnext+1] = dsys->r[k];
	senddusttonextCPU[5*indexnext+2] = dsys->th[k];
	senddusttonextCPU[5*indexnext+3] = dsys->vr[k];
	if (!ZZIntegrator)
	  senddusttonextCPU[5*indexnext+4] = dsys->vth[k];
	else
	  senddusttonextCPU[5*indexnext+4] = dsys->l[k];
	indexnext+=1;
      }
    }
  }
  DimToNext = 5*indexnext;
  DimToPrev = 5*indexprev;
  DimFromPrev = 0;
  DimFromNext = 0;
  MPI_Barrier (MPI_COMM_WORLD);
  if (CPU_Rank != CPU_Highest) {
    MPI_Isend(&DimToNext, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD, &req1);
    MPI_Wait (&req1, &fargostat);
    MPI_Irecv(&DimFromNext, 1, MPI_INT, CPU_Next, 20, MPI_COMM_WORLD, &req2);
    MPI_Wait (&req2, &fargostat);
  }
  if (CPU_Rank != 0) {
    MPI_Isend(&DimToPrev, 1, MPI_INT, CPU_Prev, 20, MPI_COMM_WORLD, &req2);
    MPI_Wait (&req2, &fargostat);
    MPI_Irecv(&DimFromPrev, 1, MPI_INT, CPU_Prev, 10, MPI_COMM_WORLD, &req1);
    MPI_Wait (&req1, &fargostat);
  }
  
  /* ========================== */
  /* Now let's send the data... */
  /* ========================== */
  MPI_Barrier (MPI_COMM_WORLD);
  /* Send to next CPU */
  if ( (CPU_Rank != CPU_Highest) && (DimToNext != 0) ) {
    MPI_Isend(senddusttonextCPU, DimToNext, MPI_DOUBLE, CPU_Next, 100, MPI_COMM_WORLD, &req3);
    MPI_Wait (&req3, &fargostat);
    if (debug == YES) {
      sprintf(name,"%dto%d.dat",CPU_Rank, CPU_Next);
      file = fopen(name, "a");
      for (i=0; i<DimToNext/5;i++) {
	fprintf(file,"%.12g\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",PhysicalTime,DimToNext,senddusttonextCPU[5*i],senddusttonextCPU[5*i+1],senddusttonextCPU[5*i+2],senddusttonextCPU[5*i+3],senddusttonextCPU[5*i+4]);
      }
      fclose(file);
    }
  }
  /* Receive from next CPU */
  if ( (CPU_Rank != CPU_Highest) && (DimFromNext != 0) ) {
    MPI_Irecv(recvdustfromnextCPU, DimFromNext, MPI_DOUBLE, CPU_Next, 200, MPI_COMM_WORLD, &req4);
    MPI_Wait (&req4, &fargostat);
    if (debug == YES) {
      sprintf(name,"%dfrom%d.dat",CPU_Rank, CPU_Next);
      file = fopen(name, "a");
      for (i=0; i<DimFromNext/5; i++) 
	fprintf(file,"%.12g\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",PhysicalTime,DimFromNext,recvdustfromnextCPU[5*i],recvdustfromnextCPU[5*i+1],recvdustfromnextCPU[5*i+2],recvdustfromnextCPU[5*i+3],recvdustfromnextCPU[5*i+4]);
      fclose(file);
    }
  }
  /* Send to previous CPU */
  if ( (CPU_Rank != 0) && (DimToPrev != 0) ) {
    MPI_Isend(senddusttoprevCPU, DimToPrev, MPI_DOUBLE, CPU_Prev, 200, MPI_COMM_WORLD, &req4);
    MPI_Wait (&req4, &fargostat);
    if (debug == YES) {
      sprintf(name,"%dto%d.dat",CPU_Rank, CPU_Prev);
      file = fopen(name, "a");
      for (i=0; i<DimToPrev/5;i++) {
	fprintf(file,"%.12g\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",PhysicalTime,DimToPrev,senddusttoprevCPU[5*i],senddusttoprevCPU[5*i+1],senddusttoprevCPU[5*i+2],senddusttoprevCPU[5*i+3],senddusttoprevCPU[5*i+4]);
      }
      fclose(file);
    }
  }
  /* Receive from previous CPU */
  if ( (CPU_Rank != 0) && (DimFromPrev != 0) ) {
    MPI_Irecv(recvdustfromprevCPU, DimFromPrev, MPI_DOUBLE, CPU_Prev, 100, MPI_COMM_WORLD, &req3);
    MPI_Wait (&req3, &fargostat);
    if (debug == YES) {
      sprintf(name,"%dfrom%d.dat",CPU_Rank, CPU_Prev);
      file = fopen(name, "a");
      for (i=0; i<DimFromPrev/5; i++) 
	fprintf(file,"%.12g\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",PhysicalTime,DimFromPrev,recvdustfromprevCPU[5*i],recvdustfromprevCPU[5*i+1],recvdustfromprevCPU[5*i+2],recvdustfromprevCPU[5*i+3],recvdustfromprevCPU[5*i+4]);
      fclose(file);
    }
  }
  /* =============================== */
  /* Now let's update particles locally */
  /* =============================== */
  if ( (CPU_Rank != 0) && (DimFromPrev != 0) ) {
    for (i=0; i<DimFromPrev/5; i++) {
      myk = (int)recvdustfromprevCPU[5*i];
      dsys->r[myk] = recvdustfromprevCPU[5*i+1];
      dsys->th[myk] = recvdustfromprevCPU[5*i+2];
      dsys->vr[myk] = recvdustfromprevCPU[5*i+3];
      if (!ZZIntegrator)
	dsys->vth[myk] = recvdustfromprevCPU[5*i+4];
      else
	dsys->l[myk] = recvdustfromprevCPU[5*i+4];
    }
  }
  if ( (CPU_Rank != CPU_Highest) && (DimFromNext != 0) ) {
    for (i=0; i<DimFromNext/5; i++) {
      myk = (int)recvdustfromnextCPU[5*i];
      dsys->r[myk] = recvdustfromnextCPU[5*i+1];
      dsys->th[myk] = recvdustfromnextCPU[5*i+2];
      dsys->vr[myk] = recvdustfromnextCPU[5*i+3];
      if (!ZZIntegrator)
	dsys->vth[myk] = recvdustfromnextCPU[5*i+4];
      else
	dsys->l[myk] = recvdustfromnextCPU[5*i+4];
    }
  }
}

void UpdateDustVelocities(dsys, Plsys, gasvr, gasvt, gasdens, dustpcdens, timestep)   
     DustSystem *dsys;
     PlanetarySystem *Plsys;
     PolarGrid *gasvr, *gasvt, *gasdens, *dustpcdens;
     real timestep;
{
  int k, i;
  real FgravR, FgravTh, d;
  real xp, yp, rp, tp, mp, rd, td, vrd, vtd, ld;
  real vrg, vtg, lg, St, stoptime;
  real PotPlan, OmegaPlan, eps, numR, numTh, den;
  int LocalDustFeelDisk;
  extern boolean IsDisk, DustFeelDisk, DustFeelPlanets, DustFeelSG, Indirect_Term, SelfGravity, ZZIntegrator, ShortFrictionTimeApproximation, DustGrowth;
  
  OmegaPlan = PotPlan = 0.0;

  /* NEW (April 2016): very simple model for dust growth due to
     Brownian motion between dust particles. All CPUs increase dust
     size so that there is no need to communicate particles size when
     particles change CPU. That's why I choose DUSTGROWTHPARAMETER to
     be constant... */
  if (DustGrowth) {
    DustMassTaper = (PhysicalTime-PhysicalTimeInitial)/(DUSTMASSTAPER*2.0*M_PI);
    DustMassTaper = (DustMassTaper > 1.0 ? 1.0 : pow(sin(DustMassTaper*M_PI/2.0),2.0));
    for (k = 0; k < NBPART; k++) {      
      dsys->dustsize[k] = dsys->dustsize_init[k] + (1e3*dsys->dustsize_init[k])*DustMassTaper;
    }
  }
  
  /* Call interpolation function to get the gas density, velocity
     vector and sound speed at particle's position */
  if (IsDisk == YES)
    interpolation(dsys, gasvr, gasvt, gasdens, dustpcdens, timestep);
  
  for (k = 0; k < NBPART; k++) {
    rd  = dsys->r[k];
    if ( (rd >= Rinf[Zero_or_active]) && (rd < Rsup[Max_or_active-1]) ) {
      td  = dsys->th[k];
      vrd = dsys->vr[k];
      vtd = dsys->vth[k];
      ld  = dsys->l[k];
      vrg = dsys->gasvratpc[k];
      vtg = dsys->gasvtatpc[k];
      St  = dsys->stokesnb[k];
      /* If the particle's leaves the grid, its velocity update no
	 longer takes gas drag into account, but only the primary's
	 gravity */
      LocalDustFeelDisk = DustFeelDisk;
      if ( (rd < RMIN) || (rd > RMAX) )
	LocalDustFeelDisk = 0;

      /* Compare particle's stopping time (ts) with hydrodynamical
	 timestep (dt): if ts < dt, update particle's velocities with
	 the short-friction time approximation of Johansen & Klahr,
	 otherwise use semi-implicite integrator below */
      stoptime = St*pow(rd,1.5);
      if ( (stoptime < timestep) && (ShortFrictionTimeApproximation) ) {
	// =================================
	// short-friction time approximation
	// CB (03/2020): expression could account for feedback?
	// =================================
	dsys->vr[k]  = vrg + stoptime*dsys->rad_gradp[k]/dsys->gasdensatpc[k];
	dsys->vth[k] = vtg + stoptime*dsys->azi_gradp[k]/dsys->gasdensatpc[k];
	dsys->l[k] = rd*dsys->vth[k];
	//printf ("SFTA: vr = %lg, vth = %lg, r = %lg, th = %lg\n",dsys->vr[k],dsys->vth[k],rd,td);
      } else {
	// =========================
	// semi-implicite integrator
	// =========================
	/* Radial and azimuthal forces from the star (direct term) */
	FgravR = -pow(rd,-2.);
	FgravTh = 0.;
      
	/* Particular case of a two-body problem between the central
	   object and the dust particle: we need this to maintain an
	   actual numerical equilibrium, otherwise the integrator will go
	   unstable */
	if ( (fabs( FgravR + vtd*vtd/rd ) < 1e-15) && (!DustFeelPlanets) && (!DustFeelDisk)  && (!DustFeelSG) )
	  FgravR = -vtd*vtd/rd;

	/* Keep track of the particles effective gravity accelerations */
	/* Radial gravity acceleration includes direct radial
	   acceleration from central star + centrifugal acceleration */
	dsys->rad_geff_acc[k] = FgravR + vtd*vtd/rd;
	/* Azimuthal gravity acceleration includes direct azimuthal
	   acceleration from central star (=0 here) - vr.vphi/r term */
	dsys->azi_geff_acc[k] = FgravTh - vrd*vtd/rd;

	/* Add the radial and azimuthal forces from the planets (direct
	   terms) */
	if (DustFeelPlanets) {
	  for (i = 0; i < Plsys->nb; i++) {
	    xp = Plsys->x[i];
	    yp = Plsys->y[i];
	    rp = sqrt( xp*xp + yp*yp);
	    tp = atan2(yp,xp);
	    mp = Plsys->mass[i];
	    /* Believe it or not, on some architectures calculating
	       the power of a number close to unity can result in an
	       extremely slow calculation (call to slow_pow
	       function). This can be the case here when rp is close
	       to one... that's why we use the trick to divide rp by 2
	       and compensate by the extra 2^1+f factor... */
	    eps = THICKNESSSMOOTHING * ASPECTRATIO * pow(0.5*rp,1.0+FLARINGINDEX) * pow(2.0,1.0+FLARINGINDEX);
	    // eps = compute_smoothing(rp);  // smoothing length at planet's orbital radius
	    // Smoothed distance between planet and dust particle
	    d = sqrt( rd*rd + rp*rp - 2.0*rd*rp*cos(td-tp) + eps*eps );
	    FgravR  +=   mp*(rp*cos(td-tp)-rd)*pow(d,-3.0);
	    FgravTh +=  -mp*rp*sin(td-tp)*pow(d,-3.0);
	    // Quantities used to compute particle's Jacobi constant if only one planet
	    PotPlan = mp/d;
	    OmegaPlan = pow(0.5*rp,-1.5) * pow(2.0,-1.5); 
	  }
	}

	/* Add indirect term from the gas disc and the planets */
	if (Indirect_Term == YES) {
	  FgravR += dsys->rad_ind_acc[k];
	  FgravTh += dsys->azi_ind_acc[k];
	}
      
	/* Add disc's gravitational potential if the disc is self-gravitating */
	if ( (IsDisk == YES) && SelfGravity && DustFeelSG) {
	  FgravR += dsys->rad_sg_acc[k];
	  FgravTh += dsys->azi_sg_acc[k];
	}
      
	if (!ZZIntegrator) {
	  /* den is used in the implicit update of dust velocities */
	  den = 1.0 + LocalDustFeelDisk*timestep/stoptime;
	  /* Implicite update particle's vtheta */
	  numTh =  (vtd-vtg) + timestep*(FgravTh - vrd*vtd/rd);
	  dsys->vth[k] = vtg + numTh/den;	
	  /* Implicite update of particle's vr. It is important in the
	     line below to keep dsys->vth[k]! */
	  numR = (vrd-vrg) + timestep*(FgravR + (dsys->vth[k]*dsys->vth[k])/rd);
	  dsys->vr[k] = vrg + numR/den;
	} else {
	  /* ZZ Integrator (Zhu et al. 2014, ApJ, 785) gives basically
	     the exact same results as the above integrator if therein
	     timestep was changed to timestep/2. Both integrators give
	     results that are in very good agreement otherwise,
	     provided that the friction time is larger than the hydro
	     timestep. ZZ integrator doesn't behave well for very
	     short friction times whereas the above integrator works
	     well and gives same results as when the short-friction
	     time approximation scheme is used */
	  /* den is used in the implicit update of dust velocities */
	  den = 1.0 + LocalDustFeelDisk*0.5*timestep/stoptime;
	  /* Implicite update particle's specific AM: l = Rxvtheta */
	  FgravTh *= rd;  // don't remove this line!
	  numTh = timestep*(FgravTh + LocalDustFeelDisk*(rd*vtg-ld)/stoptime);
	  dsys->l[k] += numTh/den;
	  dsys->vth[k] = dsys->l[k]/rd;
	  /* Implicite update of particle's vr */
	  numR = timestep*(FgravR + 0.5*pow(rd,-3.0)*(ld*ld + dsys->l[k]*dsys->l[k]) + LocalDustFeelDisk*(vrg-vrd)/stoptime);
	  dsys->vr[k] += numR/den;
	}
	if ( (dsys->vth[k] < 0.0) || (dsys->vth[k] > 10.0) || (dsys->vr[k] < -1.0) || (dsys->vr[k] > 10.0) ) {
	  printf ("Issue with update of particle %d: rd=%lg, td=%lg, vrd=%lg, vtd=%lg, vrg=%lg, vtg=%lg, Ts=%lg, Fgravth=%lg\n",k,rd,td,vrd,vtd,vrg,vtg,stoptime,FgravTh);
	}
      }
      /* Jacobi's constant. Again, it is important in the line below to
	 keep the particles updated velocities dsys->vr[k] and
	 dsys->vth[k]! */
      dsys->jacobi[k] = (0.5*(pow(dsys->vr[k],2.) + pow(dsys->vth[k],2.))) - 1./dsys->r[k] - PotPlan - OmegaPlan*rd*dsys->vth[k];
    }
  }
}

real
myrandn_rad (real mu, real sigma)
{
  real u1, u2, w, mult;
  static real x1, x2; // need to be static
  static int call = 0;
  if (call == 1)
    {
      call = !call;
      return (mu + sigma*x2);
    }
  do
    {
      u1 = 2.0*drand48()-1.0;
      u2 = 2.0*drand48()-1.0;
      w = pow(u1,2.0) + pow(u2,2.0);
    }
  while (w >= 1.0 || w == 0.0);
  mult = sqrt ((-2.0*log(w))/w);
  x1 = u1*mult;
  x2 = u2*mult;
  call = !call;
  return (mu + sigma*x1);
}

real
myrandn_phi (real mu, real sigma)
{
  real u1, u2, w, mult;
  static real x1, x2; // need to be static
  static int call = 0;
  if (call == 1)
    {
      call = !call;
      return (mu + sigma*x2);
    }
  do
    {
      u1 = 2.0*drand48()-1.0;
      u2 = 2.0*drand48()-1.0;
      w = pow(u1,2.0) + pow(u2,2.0);
    }
  while (w >= 1.0 || w == 0.0);
  mult = sqrt ((-2.0*log(w))/w);
  x1 = u1*mult;
  x2 = u2*mult;
  call = !call;
  return (mu + sigma*x1);
}
