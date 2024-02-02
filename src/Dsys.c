/* Contains allocation, initialization of dust particles, and
   interpolation function to compute gas quantities at particles
   positions */

#include "mp.h"

DustSystem *AllocDustSystem (nb)
int nb;
{
  real *r,*ds,*dsi,*the,*vr,*vth,*l;
  real *rad_geff_acc,*azi_geff_acc,*rad_drag_acc,*azi_drag_acc;
  real *rad_ind_acc,*azi_ind_acc,*rad_sg_acc,*azi_sg_acc, *rad_gradp, *azi_gradp;
  real *gasvratpc,*gasvtatpc,*gasdensatpc,*gascsatpc,*jacobi, *stokesnb;
  real *senddusttoprevCPU, *senddusttonextCPU, *recvdustfromprevCPU, *recvdustfromnextCPU;
  DustSystem *sys;
  int i;

  masterprint ("Initializing the %d dust particles\n",nb);

  sys = (DustSystem *)malloc(sizeof(DustSystem));
  if (sys == NULL) {
    fprintf(stderr, "Not enough memory to allocate dust system in Dsys.c... Aborting!\n");
    prs_exit (1);
  }
  ds = (real *)malloc(sizeof(real)*(nb+1));
  dsi = (real *)malloc(sizeof(real)*(nb+1));
  r = (real *)malloc(sizeof(real)*(nb+1));
  the = (real *)malloc(sizeof(real)*(nb+1));
  vr = (real *)malloc(sizeof(real)*(nb+1));
  vth = (real *)malloc(sizeof(real)*(nb+1));
  l = (real *)malloc(sizeof(real)*(nb+1));
  rad_gradp = (real *)malloc(sizeof(real)*(nb+1));
  azi_gradp = (real *)malloc(sizeof(real)*(nb+1));
  rad_geff_acc = (real *)malloc(sizeof(real)*(nb+1));
  azi_geff_acc = (real *)malloc(sizeof(real)*(nb+1));
  rad_drag_acc = (real *)malloc(sizeof(real)*(nb+1));
  azi_drag_acc = (real *)malloc(sizeof(real)*(nb+1));
  rad_ind_acc = (real *)malloc(sizeof(real)*(nb+1));
  azi_ind_acc = (real *)malloc(sizeof(real)*(nb+1));
  rad_sg_acc = (real *)malloc(sizeof(real)*(nb+1));
  azi_sg_acc = (real *)malloc(sizeof(real)*(nb+1));
  gasvratpc = (real *)malloc(sizeof(real)*(nb+1));
  gasvtatpc = (real *)malloc(sizeof(real)*(nb+1));
  gasdensatpc = (real *)malloc(sizeof(real)*(nb+1));
  gascsatpc = (real *)malloc(sizeof(real)*(nb+1));
  jacobi = (real *)malloc(sizeof(real)*(nb+1));
  stokesnb = (real *)malloc(sizeof(real)*(nb+1));
  senddusttoprevCPU = (real *)malloc(sizeof(real)*(5*nb));
  senddusttonextCPU = (real *)malloc(sizeof(real)*(5*nb));
  recvdustfromprevCPU = (real *)malloc(sizeof(real)*(5*nb));
  recvdustfromnextCPU = (real *)malloc(sizeof(real)*(5*nb));
  
  if ((ds==NULL) || (dsi==NULL) || (r==NULL) || (the==NULL) || (vr == NULL) || (vth==NULL) || (l==NULL) || (rad_ind_acc==NULL) || (azi_ind_acc==NULL) || (rad_gradp==NULL) || (azi_gradp==NULL) || (rad_sg_acc==NULL) || (azi_sg_acc==NULL) || (rad_geff_acc==NULL) || (azi_geff_acc==NULL) || (rad_drag_acc==NULL) || (azi_drag_acc==NULL) || (gasvratpc==NULL) ||(gasvtatpc==NULL) || (gasdensatpc==NULL) || (gascsatpc==NULL) || (jacobi == NULL) || (stokesnb == NULL) || (senddusttoprevCPU==NULL) || (senddusttonextCPU==NULL) || (recvdustfromprevCPU==NULL) || (recvdustfromnextCPU==NULL)) {
    fprintf (stderr, "Not enough memory in Dsys.c... aborting! \n");
    prs_exit (1);
  }
  
  sys->dustsize = ds;
  sys->dustsize_init = dsi;
  sys->r = r;
  sys->th = the;
  sys->vr = vr;
  sys->vth = vth;
  sys->l = l;
  sys->rad_geff_acc = rad_geff_acc;
  sys->azi_geff_acc = azi_geff_acc;
  sys->rad_gradp = rad_gradp;
  sys->azi_gradp = azi_gradp;
  sys->rad_drag_acc = rad_drag_acc;
  sys->azi_drag_acc = azi_drag_acc;
  sys->rad_ind_acc = rad_ind_acc;
  sys->azi_ind_acc = azi_ind_acc;
  sys->rad_sg_acc = rad_sg_acc;
  sys->azi_sg_acc = azi_sg_acc;
  sys->gasvratpc = gasvratpc;
  sys->gasvtatpc = gasvtatpc;
  sys->gasdensatpc = gasdensatpc;
  sys->gascsatpc = gascsatpc;
  sys->jacobi = jacobi;
  sys->stokesnb = stokesnb;
  sys->senddusttoprevCPU = senddusttoprevCPU;
  sys->senddusttonextCPU = senddusttonextCPU;
  sys->recvdustfromprevCPU = recvdustfromprevCPU;
  sys->recvdustfromnextCPU = recvdustfromnextCPU;

  for (i=0; i<nb; i++) {
    ds[i] = dsi[i] = r[i] = the[i] = vr[i] = vth[i] = l[i] = rad_ind_acc[i] = azi_ind_acc[i] = rad_geff_acc[i] = azi_geff_acc[i] = rad_gradp[i] = azi_gradp[i] = rad_drag_acc[i] = azi_drag_acc[i] = rad_sg_acc[i] = azi_sg_acc[i] = gasvratpc[i] = gasvtatpc[i] = gasdensatpc[i] = gascsatpc[i] = jacobi[i] = stokesnb[i] = senddusttoprevCPU[i] = senddusttonextCPU[i] = recvdustfromprevCPU[i] = recvdustfromnextCPU[i] = 0.0;
  }
  
  return sys;
}


void FreeDust (sys)
     DustSystem *sys;
{
  free (sys->dustsize);
  free (sys->dustsize_init);
  free (sys->r);
  free (sys->th);
  free (sys->vr);
  free (sys->vth);
  free (sys->l);
  free (sys->gasvratpc);
  free (sys->gasvtatpc);
  free (sys->rad_geff_acc);
  free (sys->azi_geff_acc);
  free (sys->rad_gradp);
  free (sys->azi_gradp);
  free (sys->rad_drag_acc);
  free (sys->azi_drag_acc);
  free (sys->rad_ind_acc);
  free (sys->azi_ind_acc);
  free (sys->rad_sg_acc);
  free (sys->azi_sg_acc);
  free (sys->gasdensatpc);
  free (sys->gascsatpc);
  free (sys->jacobi);
  free (sys->stokesnb);
  free (sys->senddusttoprevCPU);
  free (sys->senddusttonextCPU);
  free (sys->recvdustfromprevCPU);
  free (sys->recvdustfromnextCPU);
  free (sys);
}


DustSystem *InitDustSystem ()
{
  DustSystem *sys;
  int i, k=0;
  real C1, C2;
  real radius, azimuth, vrad, vtheta, ts, dsize;
  int ipl;
  real dist, ri, rip1, dr, sgacc;
  extern boolean Restart, RestartWithNewDust;
  extern int NbRestart, SelfGravity, SGZeroMode;
  FILE *input;
  char s[512], nm[512], *s1, filename[512];

  sys = AllocDustSystem(NBPART);

  /* Minimum and maximum radii between which dust particles are set
     initially */
  if (RMINDUST == 0.0) 
    RMINDUST = RMIN;
  if (RMAXDUST == 0.0) 
    RMAXDUST = RMAX;
  
  /* Check consistency between choices of SIZEMINPART, SIZEMAXPART and
     SIZEPARTSLOPE */
  if ( (SIZEMINPART == SIZEMAXPART) && (SIZEPARTSLOPE != 0.0) ) {
    mastererr ("You cannot have SIZEMINPART = SIZEMAXPART but SIZEPARTSLOPE != 0!\n");
    mastererr ("Aborted. Please check your parameters in the .par file\n");
    prs_exit(1);
  }
  
  if ( (!Restart) || (RestartWithNewDust) ) {
    /* Only CPH_Highest sorts out particle's size, radius and
       azimuth. For the particles size, we assume their probability
       distribution follows a power-law with exponent -SIZEPARTSLOPE
       between minimum size SIZEMINPART and maximum size SIZEMINPART
       (all parameters are set in the .par parameter file). For the
       particles radius, we assume their probability distribution
       follows a power-law with exponent -DUSTSLOPE which may differ
       from that of the background initial gas density. The
       prob. distribution int_Rmin^R sigma_dust(r)dr / int_Rmin^Rmax
       sigma_dust(r)dr is sorted out randomly uniformly between 0 and
       1, from which we infer the particle's R. Azimuths are sorted
       out randomly uniformly between min and max azimuths of cells
       interfaces. */
    if (CPU_Rank == CPU_Highest) {
      srand48(time(NULL));
      /* SORT OUT PARTICLES SIZE IN CODE UNITS */
      if ( SIZEMINPART != SIZEMAXPART) {
	if (SIZEPARTSLOPE != 1.0) {
	  C1 = pow(SIZEMAXPART,1.0-SIZEPARTSLOPE);
	  C2 = pow(SIZEMINPART,1.0-SIZEPARTSLOPE);
	  for (i=0; i<NBPART ;i++)
	    sys->dustsize[i] = pow(C2 + (C1-C2)*drand48(),1.0/(1.0-SIZEPARTSLOPE))/unit_length ;
	} else {
	  for (i=0; i<NBPART ;i++)
	    sys->dustsize[i] = SIZEMINPART*exp(drand48()*log(SIZEMAXPART/SIZEMINPART))/unit_length;
	}
      } else { 
	for (i=0; i<NBPART ;i++)
	  sys->dustsize[i] = SIZEMINPART/unit_length;
      }
      /* SORT OUT PARTICLES RADIUS */
      if (DUSTSLOPE != 1.0) {
	C1 = pow(RMAXDUST,1.0-DUSTSLOPE);
	C2 = pow(RMINDUST,1.0-DUSTSLOPE);
	for (i=0; i<NBPART ;i++)
	  sys->r[i] = pow(C2 + (C1-C2)*drand48(),1.0/(1.0-DUSTSLOPE));
      } else {
	for (i=0; i<NBPART ;i++)
	  sys->r[i] = RMINDUST * exp(drand48()*log(RMAXDUST/RMINDUST));
      }
      /* SORT OUT PARTICLES AZIMUTH */
      for (i=0; i<NBPART; i++)
	sys->th[i] = AziInf[0] + drand48()*(AziSup[NSEC-1]-AziInf[0]);
    }
    /* CPU_Highest communicates all particles's positions to all other CPUS */
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (sys->dustsize, NBPART, MPI_DOUBLE,CPU_Highest, MPI_COMM_WORLD);
    MPI_Bcast (sys->r, NBPART, MPI_DOUBLE,CPU_Highest, MPI_COMM_WORLD);
    MPI_Bcast (sys->th, NBPART, MPI_DOUBLE,CPU_Highest, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    /* We initialise particles with zero radial velocities and
       Keplerian azimuthal velocities. NEW: We now account for the
       self-gravitating radial acceleration of the gas in setting the
       initial azimuthal velocity of the particles. Just like for
       planets, when any. */
    if (SelfGravity) {
      if ( !SGZeroMode )
	mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	GLOBAL_AxiSGAccr = SG_Accr;
    }
    for(i=0; i<NBPART; i++) {
      sys->dustsize_init[i] = sys->dustsize[i];
      if (!SelfGravity)
	sys->vth[i] = pow(sys->r[i],-0.5);
      else {
	dist = sys->r[i];
	ipl = 0;
	while ( (GlobalRmed[ipl] <= dist) && (ipl < GLOBALNRAD-2) ) ipl++;
	ri = GlobalRmed[ipl];
	rip1 = GlobalRmed[ipl+1];
	dr = rip1 - ri;
	sgacc = (dist - ri)*GLOBAL_AxiSGAccr[ipl+1] + (rip1 - dist)*GLOBAL_AxiSGAccr[ipl];
	sgacc /= dr;
	sys->vth[i] = pow(sys->r[i],-0.5)*sqrt (1.0-dist*dist*sgacc);
      }
      sys->l[i] = sys->r[i]*sys->vth[i];
    }
  } else {
    /* Restart case: each CPU reads file dustsystatxx.dat */
    masterprint ("Restarting dust particles...\n");
    sprintf (filename, "%sdustsystat%d.dat", OUTPUTDIR, NbRestart);
    input = fopen (filename, "r");
    if (input == NULL) {
      fprintf (stderr, "WARNING ! Can't read %s to reinitialize particles. Aborting.\n", filename); 
      prs_exit(1);
    }
    while (fgets(s, NBPART, input) != NULL) {
      sscanf(s, "%lg %lg %lg %lg %lg %lg", &radius, &azimuth, &vrad, &vtheta, &ts, &dsize);
      sys->r[k] = (real)radius;
      sys->th[k] = (real)azimuth;
      sys->vr[k] = (real)vrad;
      sys->vth[k] = (real)vtheta;
      sys->stokesnb[k] = (real)ts;
      sys->dustsize[k] = (real)dsize/unit_length;
      sys->dustsize_init[k] = sys->dustsize[k];
      sys->l[k] = sys->r[k]*sys->vth[k];
      k++;
    }
    NBPART = k;
    printf ("Number of dust particles at restart is %d\n",NBPART);
  }
  return sys;
}


void interpolation(sys, gasvr, gasvt, gasdens, dustpcdens, timestep)
     DustSystem *sys;
     PolarGrid *gasvr, *gasvt, *gasdens, *dustpcdens;
     real timestep;
{
  int ip, jp, myip, myjp, i, j, myj, l, m, k, nr, ns;
  int imin, imax, jmin, jmax;
  int size_com, o, oo, one_if_odd, alloc_size;
  real rp, tp, vrp, vtp, mp, wr, wt, myazimuth;
  real dr, dphi, delta_r, delta_phi, delta, nu;
  real St, Ts;
  real Cdrag, Re, k_D, f_D, Kn, lambda, densg_cgs, Mach, dust_density_codeunits;
  real vrg, vtg, csg, densg, mg, dv2, rhop;
  real *vr, *vt, *dens, *cs;
  real *radindacc, *aziindacc;
  real *radsgacc, *azisgacc;
  real *radfbacc, *azifbacc, *fbedot;
  real *radgradp, *azigradp;
  real *dustdensity;
  extern boolean IsDisk, Indirect_Term, EnergyEquation, SelfGravity, LogGrid;
  extern boolean NGPInterpolation, CICInterpolation, TSCInterpolation, DustFeedback, DustFeelSGZeroMode, ShortFrictionTimeApproximation;
  MPI_Request req1, req2;
  static real *SendBufferPrev, *SendBufferNext;
  static real *RecvBufferPrev, *RecvBufferNext;
  static boolean allocate=YES;

  nr   = gasvr->Nrad;
  ns   = gasvr->Nsec;
  vr   = gasvr->Field;
  vt   = gasvt->Field;
  dens = gasdens->Field;
  cs   = SoundSpeed->Field; 
  radindacc = RadIndAcc->Field;
  aziindacc = AziIndAcc->Field;
  radsgacc = RadSGAcc->Field;
  azisgacc = AziSGAcc->Field;
  radfbacc = RadFBAcc->Field;
  azifbacc = AziFBAcc->Field;
  fbedot   = FBedot->Field;
  radgradp = RadGradP->Field;
  azigradp = AziGradP->Field;
  dustdensity = dustpcdens->Field;
  
  /* Convert particle's mean density from g/cm^3 to code units */
  dust_density_codeunits = RHOPART*1000.0*pow(unit_length, 3.)/unit_mass;
  
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l= j + i*ns;
      /* It is necessary to calculate the gas azimuthal velocity at the
	 particles locations in the fixed frame */
      vt[l] += Rmed[i]*OmegaFrame;
      if (DustFeedback) {
	radfbacc[l] = 0.0;
	azifbacc[l] = 0.0;
	fbedot[l]   = 0.0;
      }
      dustdensity[l] = 0.0;
    }
  }
  
  /* Grid spacing along radial and azimuthal directions */
  if (!LogGrid)
    delta_r = Rsup[0]-Rinf[0];        // arithmetic radial spacing
  else
    delta_r = log(Rsup[0]/Rinf[0]);   // logarithmic radial spacing
  delta_phi = AziSup[0]-AziInf[0];

  /* Minimum_Stopping_Tim: minimum stopping time of the
     particles. Here it is initialized to an arbitrarilly large
     number */
  Minimum_Stopping_Time = 1e6;

  /* Particular case where particles only feel the azisymmetric part
     of the disc's self-gravitating potential */
  if (DustFeelSGZeroMode)
    mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);

  /* ======================== */
  /* Loop over dust particles */
  /* ======================== */
  for (k=0; k<NBPART; k++) {
    sys->gasvratpc[k] = 0.0;
    sys->gasvtatpc[k] = 0.0;
    sys->rad_gradp[k] = 0.0;
    sys->azi_gradp[k] = 0.0;
    sys->rad_ind_acc[k] = 0.0;
    sys->rad_sg_acc[k] = 0.0;
    sys->azi_ind_acc[k] = 0.0;
    sys->azi_sg_acc[k] = 0.0;
    sys->gasdensatpc[k] = 0.0;
    sys->gascsatpc[k] = 0.0;
    rp = sys->r[k];
    /* Find out which CPU handles each particle */
    if ( (rp >= Rinf[Zero_or_active]) && (rp < Rsup[Max_or_active-1]) ) {
      /* ----------------------------------------- */
      /* 1) First find in which cell the particle is: (ip,jp) */
      /* ----------------------------------------- */
      ip = Zero_or_active;
      while ( rp >= Rinf[ip] ) ip++;
      ip--;
      tp = sys->th[k];
      jp = floor(ns*(tp-AziInf[0])/(AziSup[ns-1]-AziInf[0]));
      if (jp == ns)
	jp = 0;
      if (jp == -1) {
	printf ("jp = -1 in interpolation!");
      }
      /* --------------------------------------- */
      /* 2a) Get gas vrad at particle's location */
      /* --------------------------------------- */
      if (NGPInterpolation) {
	// Nearest-Grid-Point interpolation
	myip = ip;
	myjp = jp;
	imin = myip;
	imax = imin+1;
	jmin = myjp;
	jmax = jmin;
      }
      if (CICInterpolation) {
	// Cloud-in-Cell interpolation (bilinear)
	if (tp < Azimuth[jp])
	  myjp = jp-1;
	else
	  myjp = jp;
	myip = ip;
	imin = myip;
	imax = myip+1;
	jmin = myjp;
	jmax = myjp+1;
      }
      if (TSCInterpolation) {
	// Triangular-Shaped Cloud interpolation (quadratic splines)
	/*
	if (rp > Rmed[ip])
	  myip = ip+1;
	else 
	  myip = ip;
	*/ 
	/* We opt for a simplified radial index to avoid more comm between arrays... */
	myip = ip;
	myjp = jp;
	imin = myip-1;
	imax = myip+1;
	jmin = myjp-1;
	jmax = myjp+1;
      }
      for (j=jmin; j<=jmax; j++) {
	if ( (j >= 0) && (j < ns) ) {
	  myj = j;
	  myazimuth = Azimuth[j];
	} else {
	  if (j < 0) {
	    myj = j+ns;
	    myazimuth = Azimuth[0]+j*delta_phi;
	  }
	  if (j >= ns) {
	    myj = j-ns;
	    myazimuth = Azimuth[ns-1]+(j+1-ns)*delta_phi;
	  }
	}
	dphi = fabs(tp-myazimuth);
	for (i=imin; i<=imax; i++) {
	  if ( (i < 0) || (i > nr-1 ) ) {
	    dr = 0.0;
	  } else {
	    if (!LogGrid)
	      dr = fabs(rp-Rinf[i]);
	    else
	      dr = fabs(log(rp/Rinf[i]));
	  }
	  /* Nearest Grid Point (NGP) interpolation */
	  if (NGPInterpolation) {
	    if (dr < 0.5*delta_r) 
	      wr = 1.0;
	    else
	      wr = 0.0;
	    wt = 1.0;
	  }
	  /* Cloud-In-Cell (CIC) a.k.a. bilinear interpolation */
	  if (CICInterpolation) {
	    wr = 1.0-dr/delta_r;
	    wt = 1.0-dphi/delta_phi;
	  }
	  /* Triangular-Shaped Cloud (TSC) a.k.a. quadratic spline interpolation */
	  if (TSCInterpolation) {
	    if (dr < 0.5*delta_r) 
	      wr = 0.75-dr*dr/delta_r/delta_r;
	    if ( (dr >= 0.5*delta_r) && (dr <= 1.5*delta_r) )
	      wr = 0.5*(1.5-dr/delta_r)*(1.5-dr/delta_r);
	    if (dr > 1.5*delta_r)
	      wr = 0.0;
	    if (dphi < 0.5*delta_phi) 
	      wt = 0.75-dphi*dphi/delta_phi/delta_phi;
	    if ( (dphi >= 0.5*delta_phi) && (dphi <= 1.5*delta_phi) )
	      wt = 0.5*(1.5-dphi/delta_phi)*(1.5-dphi/delta_phi);
	    if (dphi > 1.5*delta_phi)
	      wt = 0.0;
	  }
	  if ( (wr<0.0) || (wr>1.0) || (wt<0.0) || (wt>1.0) )
	    printf ("Pb in step 2a in Dsys.c for particle %d, wr=%lg and wt=%lg; rp=%lg, i=%d, ip=%d, nr=%d, dr=%lg; tp=%lg, j=%d, jp=%d, dphi=%lg\n",k,wr,wt,rp,i,ip,nr,dr,tp,j,jp,dphi);
	  if ( (i >= 0) && (i <= nr-1) ) {
	    /* Now calculate gas radial velocity and radial pressure
	       gradient interpolated at particle's location */
	    sys->gasvratpc[k] += (wr*wt*vr[i*ns+myj]);
	    sys->rad_gradp[k] += (wr*wt*radgradp[i*ns+myj]);
	    if (Indirect_Term == YES)
	      sys->rad_ind_acc[k] += (wr*wt*radindacc[i*ns+myj]);
	    if (SelfGravity) {
	      if (!DustFeelSGZeroMode) {
		sys->rad_sg_acc[k] += (wr*wt*radsgacc[i*ns+myj]);
	      }
	      else {
		/* Particular case where particles only feel the
		   azisymmetric part of the disc's self-gravitating
		   potential */
		sys->rad_sg_acc[k] += (wr*wt*GLOBAL_AxiSGAccr[i+IMIN]);
	      }
	    }
	  }
	}
      }
      /* ----------------------------------------- */
      /* 2b) Get gas vtheta at particle's location */
      /* ----------------------------------------- */
      if (NGPInterpolation) {
	myip = ip;
	myjp = jp;
	imin = myip;
	imax = imin;
	jmin = myjp;
	jmax = jmin+1;
      }
      if (CICInterpolation) {
	if (rp < Rmed[ip])
	  myip = ip-1;
	else
	  myip = ip;
	myjp = jp;
	imin = myip;
	imax = myip+1;
	jmin = myjp;
	jmax = myjp+1;
      }
      if (TSCInterpolation) {
	if (tp > Azimuth[jp])
	  myjp = jp+1;
	else
	  myjp = jp;
	myip = ip;
	imin = myip-1;
	imax = myip+1;
	jmin = myjp-1;
	jmax = myjp+1;
      }
      for (j=jmin; j<=jmax; j++) {
	if ( (j >= 0) && (j < ns) ) {
	  myj = j;
	  myazimuth = AziInf[j];
	} else {
	  if (j < 0) {
	    myj = j+ns;
	    myazimuth = AziInf[0]+j*delta_phi;
	  }
	  if (j >= ns) {
	    myj = j-ns;
	    myazimuth = AziInf[ns-1]+(j+1-ns)*delta_phi;
	  }
	}
	dphi = fabs(tp-myazimuth);
	for (i=imin; i<=imax; i++) {	  
	  if ( (i < 0) || (i > nr-1 ) ) {
	    dr = 0.0;
	  } else {
	    if (!LogGrid)
	      dr = fabs(rp-Rmed[i]);
	    else
	      dr = fabs(log(rp/Rmed[i]));
	  }
	  /* Nearest Grid Point (NGP) interpolation */
	  if (NGPInterpolation) {
	    wr = 1.0;
	    if (dphi < 0.5*delta_phi) 
	      wt = 1.0;
	    else
	      wt = 0.0;	    
	  }
	  /* Cloud-In-Cell (CIC) a.k.a. bilinear interpolation */
	  if (CICInterpolation) {
	    wr = 1.0-dr/delta_r;
	    wt = 1.0-dphi/delta_phi;
	  }
	  /* Triangular-Shaped Cloud (TSC) a.k.a. quadratic spline interpolation */
	  if (TSCInterpolation) {
	    if (dr < 0.5*delta_r) 
	      wr = 0.75-dr*dr/delta_r/delta_r;
	    if ( (dr >= 0.5*delta_r) && (dr <= 1.5*delta_r) )
	      wr = 0.5*(1.5-dr/delta_r)*(1.5-dr/delta_r);
	    if (dr > 1.5*delta_r)
	      wr = 0.0;
	    if (dphi < 0.5*delta_phi) 
	      wt = 0.75-dphi*dphi/delta_phi/delta_phi;
	    if ( (dphi >= 0.5*delta_phi) && (dphi <= 1.5*delta_phi) )
	      wt = 0.5*(1.5-dphi/delta_phi)*(1.5-dphi/delta_phi);
	    if (dphi > 1.5*delta_phi)
	      wt = 0.0;
	  }
	  if ( (wr<0.0) || (wr>1.0) || (wt<0.0) || (wt>1.0) )
	    printf ("Pb in step 2b in Dsys.c for particle %d, wr=%lg and wt=%lg; rp=%lg, i=%d, ip=%d, nr=%d, dr=%lg, \n",k,wr,wt,rp,i,ip,nr,dr);
	  if ( (i >= 0) && (i <= nr-1) ) {
	    /* Now calculate gas vr interpolated at particle's location */
	    sys->gasvtatpc[k] += (wr*wt*vt[i*ns+myj]);
	    sys->azi_gradp[k] += (wr*wt*azigradp[i*ns+myj]);
	    if (Indirect_Term == YES)
	      sys->azi_ind_acc[k] += (wr*wt*aziindacc[i*ns+myj]);
	    if (SelfGravity) {
	      if (!DustFeelSGZeroMode) {
		sys->azi_sg_acc[k] += (wr*wt*azisgacc[i*ns+myj]);
	      }
	    }
	  }
	}
      }
      if ( (sys->gasvtatpc[k] < 0.0) || (sys->gasvtatpc[k] > 10.0) ) {
	printf ("0 - Issue with update of particle %d: vtg=%lg\n",k,sys->gasvtatpc[k]);
      }
      /* ---------------------------------------------------------- */
      /* 2c) Get gas density and sound speed at particle's location */
      /* ---------------------------------------------------------- */
      if (NGPInterpolation) {
	myip = ip;
	myjp = jp;
	imin = myip;
	imax = imin;
	jmin = myjp;
	jmax = jmin;
      }
      if (CICInterpolation) {
	if (rp < Rmed[ip])
	  myip = ip-1;
	else
	  myip = ip;
	if (tp < Azimuth[jp])
	  myjp = jp-1;
	else
	  myjp = jp;
	imin = myip;
	imax = myip+1;
	jmin = myjp;
	jmax = myjp+1;
      }
      if (TSCInterpolation) {
	myip = ip;
	myjp = jp;
	imin = myip-1;
	imax = myip+1;
	jmin = myjp-1;
	jmax = myjp+1;
      }
      for (j=jmin; j<=jmax; j++) {
	if ( (j >= 0) && (j < ns) ) {
	  myj = j;
	  myazimuth = Azimuth[j];
	} else {
	  if (j < 0) {
	    myj = j+ns;
	    myazimuth = Azimuth[0]+j*delta_phi;
	  }
	  if (j >= ns) {
	    myj = j-ns;
	    myazimuth = Azimuth[ns-1]+(j+1-ns)*delta_phi;
	  }
	}
	dphi = fabs(tp-myazimuth);
	for (i=imin; i<=imax; i++) {
	  if ( (i < 0) || (i > nr-1 ) ) {
	    dr = 0.0;
	  } else {
	    if (!LogGrid)
	      dr = fabs(rp-Rmed[i]);
	    else
	      dr = fabs(log(rp/Rmed[i]));
	  }
	  /* Nearest Grid Point (NGP) interpolation */
	  if (NGPInterpolation) {
	    wr = 1.0;
	    wt = 1.0;
	  }
	  /* Cloud-In-Cell (CIC) a.k.a. bilinear interpolation */
	  if (CICInterpolation) {
	    wr = 1.0-dr/delta_r;
	    wt = 1.0-dphi/delta_phi;
	  }
	  /* Triangular-Shaped Cloud (TSC) a.k.a. quadratic spline interpolation */
	  if (TSCInterpolation) {
	    if (dr < 0.5*delta_r) 
	      wr = 0.75-dr*dr/delta_r/delta_r;
	    if ( (dr >= 0.5*delta_r) && (dr <= 1.5*delta_r) )
	      wr = 0.5*(1.5-dr/delta_r)*(1.5-dr/delta_r);
	    if (dr > 1.5*delta_r)
	      wr = 0.0;
	    if (dphi < 0.5*delta_phi) 
	      wt = 0.75-dphi*dphi/delta_phi/delta_phi;
	    if ( (dphi >= 0.5*delta_phi) && (dphi <= 1.5*delta_phi) )
	      wt = 0.5*(1.5-dphi/delta_phi)*(1.5-dphi/delta_phi);
	    if (dphi > 1.5*delta_phi)
	      wt = 0.0;
	    // Mignone+ 2019 weights for cylindrical radius
	    // dr corresponds to delta in Mignone+ notations
	    /*
	    if ( (i >= 0) && (i <= nr-1 ) ) {
	      if (!LogGrid) {
		delta = (rp-Rmed[i])/delta_r;
		nu = Rmed[i]/delta_r;
	      }
	      else {
		delta = (log(rp/Rmed[i]))/delta_r;
		nu = log(Rmed[i]/delta_r);
	      }
	      if (i == myip) // Wi
		wr = (delta + 3.0*nu)/(3.0*delta + 3.0*nu)*(0.75-delta*delta);
	      if (i == imin) // Wi-1
		wr = (delta + 3.0*nu - 2.0)/(3.0*delta + 3.0*nu)*0.5*(0.5-delta)*(0.5-delta);
	      if (i == imax) // Wi+1
		wr = (delta + 3.0*nu + 2.0)/(3.0*delta + 3.0*nu)*0.5*(0.5+delta)*(0.5+delta);
	      dphi = tp-myazimuth;
	      if (j == myjp) // Wj
		wt = (0.75-dphi*dphi/delta_phi/delta_phi);
	      if (j == jmax) // Wj-1
		wt = 0.5*(0.5-dphi/delta_phi)*(0.5-dphi/delta_phi);
	      if (j == jmin) // Wj+1
		wt = 0.5*(0.5+dphi/delta_phi)*(0.5+dphi/delta_phi);
	    }
	    */
	  }
	  if ( (wr<0.0) || (wr>1.0) || (wt<0.0) || (wt>1.0) )
	    printf ("Pb in step 2c in Dsys.c for particle %d, wr=%lg, wt=%lg, i=%d, j=%d\n",k,wr,wt,i,myj);
	  if ( (i >= 0) && (i <= nr-1) ) {
	    /* Now calculate gas surface density and sound speed
	       interpolated at particle's location */
	    sys->gasdensatpc[k] += (wr*wt*dens[i*ns+myj]);
	    sys->gascsatpc[k]   += (wr*wt*cs[i*ns+myj]);
	    dustdensity[i*ns+myj] += wr*wt*Particles_Mass/Surf[i];
	  }
	}
      }
      /* ------------------------------------------ */
      /* 2d) Now calculate particle's stokes number */
      /* ------------------------------------------ */
      // Defined as in Paardekooper 07 as tau_friction = Stokes nb / Omega_Kep
      // assuming rho_gas = sigma_gas / sqrt(2pi) / H, 
      // we get tau_friction = pi/2 . Cdrag . (s x rho_pc / sigma_gas)
      // with s = particle size, rho_pc = particle internal density
      // Cdrag = (3Kn + 1)^2 / (9Kn^2 f_D + 3Kn k_D)
      // Kn = lambda / 2s the Knudsen number
      // lambda is the molecular mean-free path
      // f_D = sqrt(1 + 9piMa^2 / 128), Ma the relative Mach number
      // and for k_D see Paardekooper 07!
      if (IsDisk == YES) {
	vrp = sys->vr[k];
	vtp = sys->vth[k];
	vrg = sys->gasvratpc[k];
	vtg = sys->gasvtatpc[k];
	csg = sys->gascsatpc[k];
	densg = sys->gasdensatpc[k];
	/* Molecular mean-free path ~ m_H2 H / Sigma_gas d^2 with m_H2
	   mass of H2, H pressure scale height and d ~ typical
	   diameter of a particle */
	densg_cgs = 0.1*densg*unit_mass*pow(unit_length,-2.0);  // in g cm^-2
	lambda = 3.34e-8 * pow(densg_cgs,-1.) * AspectRatio(rp) * rp; // in code units
	/* Knudsen number, note that sys->dustsize[k] is already in code units */
	Kn = 0.5*lambda/sys->dustsize[k];
	/* Gas relative Mach number at particle's position: |dV| / cs */
	Mach = sqrt( (vrp-vrg)*(vrp-vrg) + (vtp-vtg)*(vtp-vtg) ) / csg;
	/* f_D coefficient to link the subsonic and the supersonic regimes */
	f_D = sqrt( 1.0 + 9.0*M_PI*Mach*Mach/128.0 );
	/* Reynolds number */
	Re = 3.0*sqrt(PI/8.0)*Mach/Kn;
	/* k_D coefficient that describes the Stokes regime */
	if (Re <= 500.0)
	  k_D = 1.0+0.15*pow(Re,0.687);
	else {
	  if ( (Re > 500.0) && (Re <= 1500.0) )
	    k_D = 3.96e-6*pow(Re,2.4);
	  if (Re > 1500.0)
	    k_D = 0.11*Re;
	}
	/* Final expression for drag coefficient Cdrag as in Paardekooper 07 */
	Cdrag = pow(3.0*Kn+1.0,2.0) * pow(9.0*Kn*Kn*f_D + 3.0*Kn*k_D,-1.0);
	/* Dimensionless stopping time at the particle = stokes number */
	St = 0.5*M_PI*Cdrag*sys->dustsize[k]*dust_density_codeunits/densg;  // Stokes number (St)
	Ts = St*pow(sys->r[k],1.5);                 // Stopping time (Ts)
	if (Ts < Minimum_Stopping_Time)  // record the minimum stopping time
	  Minimum_Stopping_Time = Ts;
      } else {
	St = 1.0;
	Ts = 1.0;
      }
      sys->stokesnb[k] = St;
      if (St < 1e-10) {
	printf ("Pb in calculating Stokes number for particle %d: St=%lg, rp=%lg, tp=%lg, densg=%lg, vrp=%lg, vtp=%lg, vrg=%lg, vtg=%lg, csg=%lg\n",k,St,rp,tp,densg,vrp,vtp,vrg,vtg,csg);
	printf ("Suite: wr=%lg, wt=%lg, i=%d, ns=%d, myj=%d, dens=%lg",wr,wt,i,ns,myj,dens[i*ns+myj]);
      }
      
      /* ----------------------------------- */
      /* 3) Particle's gas drag acceleration */
      /* ----------------------------------- */
      if (IsDisk == YES) {
	 if ( (Ts < timestep) && (ShortFrictionTimeApproximation) ) {
	   /* short-friction time approximation for dust feedback */
	   vrp = vrg + Ts*sys->rad_gradp[k]/densg;
	   vtp = vtg + Ts*sys->azi_gradp[k]/densg;
	 } 
	 sys->rad_drag_acc[k] = -(vrp-vrg)/Ts;
	 sys->azi_drag_acc[k] = -(vtp-vtg)/Ts;
      }
      if (DustFeedback) {
	/* ------------------------------------------------------------- */
	/* 4a) Infer radial acceleration due to dust feedback on the gas */
	/* ------------------------------------------------------------- */
	if (NGPInterpolation) {
	  myip = ip;
	  myjp = jp;
	  imin = myip;
	  imax = imin+1;
	  jmin = myjp;
	  jmax = jmin;
	}
	if (CICInterpolation) {
	  if (tp < Azimuth[jp])
	    myjp = jp-1;
	  else
	    myjp = jp;
	  myip = ip;
	  imin = myip;
	  imax = myip+1;
	  jmin = myjp;
	  jmax = myjp+1;
	}
	if (TSCInterpolation) {
	  /*
	  if (rp > Rmed[ip])
	    myip = ip+1;
	  else 
	    myip = ip;
	  */
	  /* We opt for a simplified radial index to avoid more comm between arrays... */
	  myip = ip;
	  myjp = jp;
	  imin = myip-1;
	  imax = myip+1;
	  jmin = myjp-1;
	  jmax = myjp+1;
	}
	for (j=jmin; j<=jmax; j++) {
	  if ( (j >= 0) && (j < ns) ) {
	    myj = j;
	    myazimuth = Azimuth[j];
	  } else {
	    if (j < 0) {
	      myj = j+ns;
	      myazimuth = Azimuth[0]+j*delta_phi;
	    }
	    if (j >= ns) {
	      myj = j-ns;
	      myazimuth = Azimuth[ns-1]+(j+1-ns)*delta_phi;
	    }
	  }
	  dphi = fabs(tp-myazimuth);
	  for (i=imin; i<=imax; i++) {
	    if ( (i < 0) || (i > nr-1 ) ) {
	      dr = 0.0;
	    } else {
	      if (!LogGrid)
		dr = fabs(rp-Rinf[i]);
	      else
		dr = fabs(log(rp/Rinf[i]));
	    }
	    /* Nearest Grid Point (NGP) interpolation */
	    if (NGPInterpolation) {
	      if (dr < 0.5*delta_r) 
		wr = 1.0;
	      else
		wr = 0.0;
	      wt = 1.0;
	    }
	    /* Cloud-In-Cell (CIC) a.k.a. bilinear interpolation */
	    if (CICInterpolation) {
	      wr = 1.0-dr/delta_r;
	      wt = 1.0-dphi/delta_phi;
	    }
	    /* Triangular-Shaped Cloud (TSC) a.k.a. quadratic spline interpolation */
	    if (TSCInterpolation) {
	      if (dr < 0.5*delta_r) 
		wr = 0.75-dr*dr/delta_r/delta_r;
	      if ( (dr >= 0.5*delta_r) && (dr <= 1.5*delta_r) )
		wr = 0.5*(1.5-dr/delta_r)*(1.5-dr/delta_r);
	      if (dr > 1.5*delta_r)
		wr = 0.0;
	      if (dphi < 0.5*delta_phi) 
		wt = 0.75-dphi*dphi/delta_phi/delta_phi;
	      if ( (dphi >= 0.5*delta_phi) && (dphi <= 1.5*delta_phi) )
		wt = 0.5*(1.5-dphi/delta_phi)*(1.5-dphi/delta_phi);
	      if (dphi > 1.5*delta_phi)
		wt = 0.0;
	    }
	    if ( (wr < 0.0) || (wr>1.0) || (wt<0.0) || (wt>1.0) )
	      printf ("Pb in step 4a in Dsys.c for particle %d, wr=%lg and wt=%lg\n",k,wr,wt);
	    if ( (i >= 0) && (i <= nr-1) ) {
	      /* Now calculate feedback radial acceleration */
	      mp = Particles_Mass;                     // (super-)particle's mass
	      mg = dens[i*ns+myj]*Surf[i];             // gas mass in cell
	      radfbacc[i*ns+myj] -= (wr*wt*sys->rad_drag_acc[k])*mp/mg;
	    }
	  }
	}
	/* ---------------------------------------------------------------- */
	/* 4b) Infer azimuthal acceleration due to dust feedback on the gas */
	/* ---------------------------------------------------------------- */
	if (NGPInterpolation) {
	  myip = ip;
	  myjp = jp;
	  imin = myip;
	  imax = imin;
	  jmin = myjp;
	  jmax = jmin+1;
	}
	if (CICInterpolation) {
	  if (rp < Rmed[ip])
	    myip = ip-1;
	  else
	    myip = ip;
	  myjp = jp;
	  imin = myip;
	  imax = myip+1;
	  jmin = myjp;
	  jmax = myjp+1;
	}
	if (TSCInterpolation) {
	  if (tp > Azimuth[jp])
	    myjp = jp+1;
	  else
	    myjp = jp;
	  myip = ip;
	  imin = myip-1;
	  imax = myip+1;
	  jmin = myjp-1;
	  jmax = myjp+1;
	}
	for (j=jmin; j<=jmax; j++) {
	  if ( (j >= 0) && (j < ns) ) {
	    myj = j;
	    myazimuth = AziInf[j];
	  } else {
	    if (j < 0) {
	      myj = j+ns;
	      myazimuth = AziInf[0]+j*delta_phi;
	    }
	    if (j >= ns) {
	      myj = j-ns;
	      myazimuth = AziInf[ns-1]+(j+1-ns)*delta_phi;
	    }
	  }
	  dphi = fabs(tp-myazimuth);
	  for (i=imin; i<=imax; i++) {
	    if ( (i < 0) || (i > nr-1 ) ) {
	      dr = 0.0;
	    } else {
	      if (!LogGrid)
		dr = fabs(rp-Rmed[i]);
	      else
		dr = fabs(log(rp/Rmed[i]));
	    }
	    /* Nearest Grid Point (NGP) interpolation */
	    if (NGPInterpolation) {
	      wr = 1.0;
	      if (dphi < 0.5*delta_phi) 
		wt = 1.0;
	      else
		wt = 0.0;
	    }
	    /* Cloud-In-Cell (CIC) a.k.a. bilinear interpolation */
	    if (CICInterpolation) {
	      wr = 1.0-dr/delta_r;
	      wt = 1.0-dphi/delta_phi;
	    }
	    /* Triangular-Shaped Cloud (TSC) a.k.a. quadratic spline interpolation */
	    if (TSCInterpolation) {
	      if (dr < 0.5*delta_r) 
		wr = 0.75-dr*dr/delta_r/delta_r;
	      if ( (dr >= 0.5*delta_r) && (dr <= 1.5*delta_r) )
		wr = 0.5*(1.5-dr/delta_r)*(1.5-dr/delta_r);
	      if (dr > 1.5*delta_r)
		wr = 0.0;
	      if (dphi < 0.5*delta_phi) 
		wt = 0.75-dphi*dphi/delta_phi/delta_phi;
	      if ( (dphi >= 0.5*delta_phi) && (dphi <= 1.5*delta_phi) )
		wt = 0.5*(1.5-dphi/delta_phi)*(1.5-dphi/delta_phi);
	      if (dphi > 1.5*delta_phi)
		wt = 0.0;
	    }
	    if ( (wr < 0.0) || (wr>1.0) || (wt<0.0) || (wt>1.0) )
	      printf ("Pb in step 4b in Dsys.c for particle %d, wr=%lg and wt=%lg\n",k,wr,wt);
	    if ( (i >= 0) && (i <= nr-1) ) {
	      /* Now calculate feedback azimuthal acceleration */
	      mp = Particles_Mass; // (super-)particle's mass
	      mg = dens[i*ns+myj]*Surf[i];             // gas mass in cell
	      azifbacc[i*ns+myj] -= (wr*wt*sys->azi_drag_acc[k])*mp/mg;
	    }
	  }
	}
	/* ------------------------------------------------------------------------------ */
	/* 4c) Infer thermal energy density increase rate due to dust feedback on the gas */
	/* ------------------------------------------------------------------------------ */
	if (EnergyEquation) {
	  if (NGPInterpolation) {
	    myip = ip;
	    myjp = jp;
	    imin = myip;
	    imax = imin;
	    jmin = myjp;
	    jmax = jmin;
	  }
	  if (CICInterpolation) {
	    if (rp < Rmed[ip])
	      myip = ip-1;
	    else
	      myip = ip;
	    if (tp < Azimuth[jp])
	      myjp = jp-1;
	    else
	      myjp = jp;
	    imin = myip;
	    imax = myip+1;
	    jmin = myjp;
	    jmax = myjp+1;
	  }
	  if (TSCInterpolation) {
	    myip = ip;
	    myjp = jp;
	    imin = myip-1;
	    imax = myip+1;
	    jmin = myjp-1;
	    jmax = myjp+1;
	  }
	  for (j=jmin; j<=jmax; j++) {
	    if ( (j >= 0) && (j < ns) ) {
	      myj = j;
	      myazimuth = Azimuth[j];
	    } else {
	      if (j < 0) {
		myj = j+ns;
		myazimuth = Azimuth[0]+j*delta_phi;
	      }
	      if (j >= ns) {
		myj = j-ns;
		myazimuth = Azimuth[ns-1]+(j+1-ns)*delta_phi;
	      }
	    }
	    dphi = fabs(tp-myazimuth);
	    for (i=imin; i<=imax; i++) {
	      if ( (i < 0) || (i > nr-1 ) ) {
		dr = 0.0;
	      } else {
		if (!LogGrid)
		  dr = fabs(rp-Rmed[i]);
		else
		  dr = fabs(log(rp/Rmed[i]));
	      }
	      /* Nearest Grid Point (NGP) interpolation */
	      if (NGPInterpolation) {
		wr = 1.0;
		wt = 1.0;
	      }
	      /* Cloud-In-Cell (CIC) a.k.a. bilinear interpolation */
	      if (CICInterpolation) {
		wr = 1.0-dr/delta_r;
		wt = 1.0-dphi/delta_phi;
	      }
	      /* Triangular-Shaped Cloud (TSC) a.k.a. quadratic spline interpolation */
	      if (TSCInterpolation) {
		if (dr < 0.5*delta_r) 
		  wr = 0.75-dr*dr/delta_r/delta_r;
		if ( (dr >= 0.5*delta_r) && (dr <= 1.5*delta_r) )
		  wr = 0.5*(1.5-dr/delta_r)*(1.5-dr/delta_r);
		if (dr > 1.5*delta_r)
		  wr = 0.0;
		if (dphi < 0.5*delta_phi) 
		  wt = 0.75-dphi*dphi/delta_phi/delta_phi;
		if ( (dphi >= 0.5*delta_phi) && (dphi <= 1.5*delta_phi) )
		  wt = 0.5*(1.5-dphi/delta_phi)*(1.5-dphi/delta_phi);
		if (dphi > 1.5*delta_phi)
		  wt = 0.0;
	      }
	      if ( (wr<0.0) || (wr>1.0) || (wt<0.0) || (wt>1.0) )
		printf ("Pb in step 4c in Dsys.c for particle %d, wr=%lg and wt=%lg\n",k,wr,wt);
	      if ( (i >= 0) && (i <= nr-1) ) {
		/* Now calculate feedback for edot gas */
		mp   = Particles_Mass; // (super-)particle's mass
		rhop = mp/Surf[i];     // dust "surface density"
		dv2 = (vrp-vrg)*(vrp-vrg) + (vtp-vtg)*(vtp-vtg);
		fbedot[i*ns+myj] += rhop*wr*wt*dv2/Ts;
	      }
	    }
	  }
	}  // end energy equation condition for dust feedback
      } // end dust feedback if condition
    } // end if condition on which CPU handles particle
  } // end loop on particles
  /* ----------------------------------------- */
  /* We put the gas azimuthal velocities back to their values in the
     corotating frame */
  /* ----------------------------------------- */
  for (i=0; i<nr; i++) {
    for (j=0; j<NSEC; j++) {
      l= j + i*NSEC;
      vt[l] -= Rmed[i]*OmegaFrame;
    }
  }
  if ( (PhysicalTime == 0) && (Minimum_Stopping_Time < 1e6) )
    printf ("At t=0, the particles' minimum stopping time is %lg in CPU %d\n", Minimum_Stopping_Time, CPU_Rank);


  /* ----------------------------------------- */
  /* We now sum the radial and azimuthal feedback accelerations
     between neighboring overlapping rings (only one ring actually) */
  /* ----------------------------------------- */

  if (CPU_Number > 1) {

    if (!DustFeedback)
      alloc_size = 1;
    else
      alloc_size = 3;
    
    if (allocate) {
      allocate = NO;
      SendBufferPrev = (real *)malloc(alloc_size*NSEC*sizeof(real));
      SendBufferNext = (real *)malloc(alloc_size*NSEC*sizeof(real));
      RecvBufferPrev = (real *)malloc(alloc_size*NSEC*sizeof(real));
      RecvBufferNext = (real *)malloc(alloc_size*NSEC*sizeof(real));
    }

    for (j=0; j<NSEC; j++) {
      SendBufferPrev[j] = dustdensity[(Zero_or_active-1)*NSEC+j];
      SendBufferNext[j] = dustdensity[Max_or_active*NSEC+j];
      if (DustFeedback) {
	SendBufferPrev[j+NSEC]   = radfbacc[(Zero_or_active-1)*NSEC+j];
	SendBufferNext[j+NSEC]   = radfbacc[Max_or_active*NSEC+j];
	SendBufferPrev[j+2*NSEC] = azifbacc[(Zero_or_active-1)*NSEC+j];
	SendBufferNext[j+2*NSEC] = azifbacc[Max_or_active*NSEC+j];
      }
    }

    MPI_Barrier (MPI_COMM_WORLD);
    
    if (CPU_Rank != CPU_Highest) {
      MPI_Isend(&SendBufferNext[0], alloc_size*NSEC, MPI_DOUBLE, CPU_Next, 10, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      MPI_Irecv(&RecvBufferNext[0], alloc_size*NSEC, MPI_DOUBLE, CPU_Next, 20, MPI_COMM_WORLD, &req2);
      MPI_Wait (&req2, &fargostat);
      for (j=0; j<NSEC; j++) {
	dustdensity[(Max_or_active-1)*NSEC+j] += RecvBufferNext[j];
	if (DustFeedback) {
	  radfbacc[(Max_or_active-1)*NSEC+j] += RecvBufferNext[j+NSEC];
	  azifbacc[(Max_or_active-1)*NSEC+j] += RecvBufferNext[j+2*NSEC];
	}
      }
    }
    if (CPU_Rank != 0) {
      MPI_Isend(&SendBufferPrev[0], alloc_size*NSEC, MPI_DOUBLE, CPU_Prev, 20, MPI_COMM_WORLD, &req2);
      MPI_Wait (&req2, &fargostat);
      MPI_Irecv(&RecvBufferPrev[0], alloc_size*NSEC, MPI_DOUBLE, CPU_Prev, 10, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (j=0; j<NSEC; j++) {
	dustdensity[Zero_or_active*NSEC+j] += RecvBufferPrev[j];
	if (DustFeedback) {
	  radfbacc[Zero_or_active*NSEC+j] += RecvBufferPrev[j+NSEC];
	  azifbacc[Zero_or_active*NSEC+j] += RecvBufferPrev[j+2*NSEC];
	}
      }
    }
  }
}



void RotateDsys (dsys, angle)	/* Rotate Dust system by angle '-angle' */
     DustSystem *dsys;
     real angle;
{
  int k;
  for (k = 0; k < NBPART; k++) {
    dsys->th[k] = dsys->th[k] - angle;
    /* Retain particle's azimuths between AziInf[0] and
       AziSup[ns-1], the first and last azimuthal interfaces of the
       grid cells. With a 2PI azimuthal extent, AziInf[0] = 0 -
       pi/ns, and AziSup[ns-1] = 2pi-pi/ns. */
    if (dsys->th[k] < AziInf[0])
      dsys->th[k] += (AziSup[NSEC-1]-AziInf[0]);  // add 2pi
    if (dsys->th[k] > AziSup[NSEC-1])
      dsys->th[k] -= (AziSup[NSEC-1]-AziInf[0]);  // subtract 2pi
  }
}
