/** \file Force.c

Contains the function used to evaluate the %force due to the disk, and
the function that writes the 'tqwk' log files. We calculate the force
on the planet by smoothly excluding a certain amount of the planet's
Hill Radius R_H. The exclusion function is equal to 0 below the
exclusion distance, equal to 1 beyond R_H, and rises as a sinus
squared in between. In absence of self-gravity, 11 exclusion distances
are considered (dimfxy is set by default to 11 in main.c), that is we
exclude 0/10, 1/10, 2/10,..., 10/10 times R_H. Each exclusion factor
is related to a unique file, e.g. tqwk0_4.dat for an exclusion of 4/10
R_H. Although the planet mass is given as an argument to
ComputeForce(), this mass is used only to specify the distance cutoff
in the case of the Hill sphere avoidance.  The force returned is a
specific force. It has therefore the dimension of an acceleration
(LT^-2).

October 2014: we no longer assign dimfxy to 11. It is now set to 2 to 
save computing time (15 to 20%). Only the case 0 and 1 times R_H are 
accounted for.
*/

#include "mp.h"

extern boolean OpenInner, NonReflecting, BM08;
extern Pair DiskOnPrimaryAcceleration;

Force *AllocateForce ()
{
  int i;
  Force *force;
  real *globalforce;
  force  = (Force *) prs_malloc (sizeof(Force));
  /* dimfxy is a global integer defined in global.h. It is set to 2
     (oct. 2014) */
  globalforce = (real *) prs_malloc (sizeof(real) * 4 * dimfxy);
  for (i = 0; i < 4*dimfxy; i++)
    globalforce[i] = 0.;
  force->GlobalForce = globalforce;
  // NEW (March 2019)
  if (BM08)
    Residual_Rho = CreatePolarGrid(NRAD, NSEC, "RRho");
  return force;
}

void FreeForce (force)
     Force *force;
{
  free (force->GlobalForce);
}


// NEW (March 2019)
void ComputeResidualDensity (Rho)
     PolarGrid *Rho;
{
  int i, j, l, ns;
  real *dens, *residual_dens, buf;
  dens = Rho->Field;
  residual_dens = Residual_Rho->Field;
  ns = Rho->Nsec;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    buf = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      buf += dens[l];
    }
    buf /= (real)ns;  // azimuthally-averaged surface density at ring i
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      residual_dens[l] = dens[l] - buf;  // residual surface density
    }
  }
}


void ComputeForce (force, Rho, x, y, rsmoothing, mass, sys)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
     PlanetarySystem *sys;
{
  int i, j, l, ns, k;
  real xc, yc, cellmass, dx, dy, distance, d2, dist2, rh, a;
  int dimexclude;
  real x0, x1, y0, y1, m0, m1, xb, yb;
  real planet_distance, cutoff;
  real InvDist3, hill_cut, hillcutfactor;
  real *fxi, *fxo, *fyi, *fyo;
  real *localforce, *globalforce;
  real *dens, *abs, *ord;
  real cutoffdist, raddistfrompla;
  /* dimfxy is a global integer defined in global.h. It is set to 2
     (oct. 2014) */
  fxi = (real *) prs_malloc (sizeof(real) * dimfxy);
  fxo = (real *) prs_malloc (sizeof(real) * dimfxy);
  fyi = (real *) prs_malloc (sizeof(real) * dimfxy);
  fyo = (real *) prs_malloc (sizeof(real) * dimfxy);
  localforce = (real *) prs_malloc (sizeof(real) * 4 * dimfxy);
  globalforce = force->GlobalForce;
  ns = Rho->Nsec;
  /* The trick below amounts to subtracting the azimuthally averaged
     density prior to the torque evaluation. This has no impact on the
     torque, but has on the angular speed of the planet and is
     required for a proper location of resonances in a non
     self-gravitating disk. See Baruteau & Masset 2008, ApJ, 678, 483
     (arXiv:0801.4413) for details. Inspired from FARGO3D source
     files! */
  if (BM08) {
    ComputeResidualDensity (Rho);
    dens = Residual_Rho->Field;
  } else {
    dens = Rho->Field;
  }
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  a = sqrt(x*x+y*y);  // star-planet distance
  rh = pow(mass/3., 1./3.)*a+1e-15;
  for ( k = 0; k < dimfxy; k++ ) {
    fxi[k] = 0.;
    fxo[k] = 0.;
    fyi[k] = 0.;
    fyo[k] = 0.;
  }
  for ( k = 0; k < 4*dimfxy; k++ ) {
    localforce[k] = 0.;
    globalforce[k] = 0.;
  }
  if (FakeSequential && (CPU_Rank > 0)) {
    MPI_Recv (&globalforce[0], 4*dimfxy, MPI_DOUBLE, CPU_Rank-1, 27, MPI_COMM_WORLD, &fargostat);
    for ( k = 0; k < dimfxy; k++ ) {
      fxi[k] = globalforce [k];
      fxo[k] = globalforce [k + dimfxy];
      fyi[k] = globalforce [k + 2*dimfxy];
      fyo[k] = globalforce [k + 3*dimfxy];
    }
  }
  if (sys->Binary[0] == YES) {
    /* Case of a binary-star system. Exclusion is done with
       respect to the barycenter of the two stars, the position of
       which is determined below */
    x0 = sys->x[0];
    x1 = sys->x[1];
    y0 = sys->y[0];
    y1 = sys->y[1];
    m0 = sys->mass[0];
    m1 = sys->mass[1];
    xb = (m0*x0 + m1*x1) / (m0+m1);
    yb = (m0*y0 + m1*y1) / (m0+m1);
  }
  
#pragma omp parallel for private(j,hill_cut,cellmass,l,xc,yc,dist2,distance,InvDist3,dx,dy) shared(fxi,fyi,fxhi,fyhi,fxo,fyo,fxho,fyho)
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      xc = abs[l];
      yc = ord[l];
      cellmass = Surf[i]*dens[l];
      dx = xc-x;
      dy = yc-y;
      d2 = dx*dx+dy*dy;
      planet_distance = sqrt(d2);
      // here rsmoothing should be substituted by your eps_p(i,j) expression!
      dist2 = d2 + rsmoothing*rsmoothing;
      distance = sqrt(dist2);
      InvDist3 = 1.0/dist2/distance;
      if (sys->Binary[0] == YES) {
	/* Remember if a binary-star is assumed, exclusion is
	   effective from the center-of-mass of the binary-star
	   system, instead of the location of each star:
	   planet_distance holds here the distance from a cell to the
	   barycenter of the two satellites */
	d2 = (xc-xb)*(xc-xb) + (yc-yb)*(yc-yb);
	planet_distance = sqrt(d2);
      }
      /* --------------- */
      for ( k = 0; k < dimfxy; k++ ) {
	if (dimfxy != 2)
	  hillcutfactor = (real)k / (real)(dimfxy-1);
	else
	  hillcutfactor = EXCLUDEHILLFACTOR;
	if ( k != 0 ) {
	  cutoff = hillcutfactor * rh;
	  /* New default exclusion function */
	  if (planet_distance/cutoff < 0.5)
	    hill_cut = 0.0;
	  else {
	    if (planet_distance > cutoff) 
	      hill_cut = 1.0;
	    else
	      hill_cut = pow(sin((planet_distance/cutoff-.5)*M_PI),2.);
	  }
	  /* Old default exclusion function */
	  //hill_cut = 1.-exp(-d2/(cutoff*cutoff));
	}
	else
	  hill_cut = 1.; // if k=0
	if (Rmed[i] < a) {
#pragma omp atomic
	  fxi[k] += G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	  fyi[k] += G*cellmass*dy*InvDist3*hill_cut;
	} else {
#pragma omp atomic
	  fxo[k] += G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	  fyo[k] += G*cellmass*dy*InvDist3*hill_cut;
	}
      }
    }
  }
  if (FakeSequential) {
    for ( k = 0; k < dimfxy; k++ ) {
      globalforce [k]            = fxi[k];
      globalforce [k + dimfxy]   = fxo[k];
      globalforce [k + 2*dimfxy] = fyi[k];
      globalforce [k + 3*dimfxy] = fyo[k];
    }
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&globalforce[0], 4*dimfxy, MPI_DOUBLE, CPU_Rank+1, 27, MPI_COMM_WORLD);
  } else {
    for ( k = 0; k < dimfxy; k++ ) {
      localforce [k]            = fxi[k];
      localforce [k + dimfxy]   = fxo[k];
      localforce [k + 2*dimfxy] = fyi[k];
      localforce [k + 3*dimfxy] = fyo[k];
    }
    MPI_Allreduce (localforce, globalforce, 4*dimfxy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential)
    MPI_Bcast (globalforce, 4*dimfxy, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  if (dimfxy != 2)
    dimexclude = (int)(EXCLUDEHILLFACTOR*(dimfxy-1));
  else
    dimexclude = 1;
  force->fx_inner    = globalforce[0];
  force->fx_ex_inner = globalforce[dimexclude];
  force->fx_outer    = globalforce[dimfxy];
  force->fx_ex_outer = globalforce[dimfxy+dimexclude];
  force->fy_inner    = globalforce[2*dimfxy];
  force->fy_ex_inner = globalforce[2*dimfxy+dimexclude];
  force->fy_outer    = globalforce[3*dimfxy];
  force->fy_ex_outer = globalforce[3*dimfxy+dimexclude];
  force->GlobalForce = globalforce;
  free (localforce);
  free (fxi);
  free (fxo);
  free (fyi);
  free (fyo);
}

real compute_smoothing (r)
     real r;
{
  real smooth;
  // disc aspect ratio scales as r^flaring_index
  // H(r) = r x h(r) 
  smooth = THICKNESSSMOOTHING * AspectRatio(r) * pow(r, 1.0+FLARINGINDEX);
  return smooth;
}


void UpdateLog (fc, psys, Rho, outputnb, time)
     Force *fc;
     PlanetarySystem *psys;
     PolarGrid *Rho;
     int outputnb;
     real time;
{
  int i, nb;
  real x, y, r, m, vx, vy, smoothing;
  real *globalforce;
  FILE *out;
  char filename[MAX1D];
  nb = psys->nb;
  
  for (i = 0; i < nb; i++) {
    x = psys->x[i];
    y = psys->y[i];
    vx = psys->vx[i];
    vy = psys->vy[i];
    r = sqrt(x*x+y*y);
    m = psys->mass[i];
    if (RocheSmoothing)
      smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
    else
      smoothing = compute_smoothing(r);
    /* dimfxy is a global integer defined in global.h. It is set to 2
       (oct. 2014) */
    ComputeForce (fc, Rho, x, y, smoothing, m, psys);
    globalforce = fc->GlobalForce;
    if (CPU_Rank == CPU_Number-1) {
      sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
      out = fopen (filename, "a");
      if (out == NULL) {
	fprintf (stderr, "Can't open %s\n", filename);
	fprintf (stderr, "Aborted.\n");
      }
      fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb, \
	       x*fc->fy_inner-y*fc->fx_inner,				\
	       x*fc->fy_outer-y*fc->fx_outer,				\
	       x*fc->fy_ex_inner-y*fc->fx_ex_inner,			\
	       x*fc->fy_ex_outer-y*fc->fx_ex_outer,			\
	       vx*fc->fx_inner+vy*fc->fy_inner,				\
	       vx*fc->fx_outer+vy*fc->fy_outer,				\
	       vx*fc->fx_ex_inner+vy*fc->fy_ex_inner,			\
	       vx*fc->fx_ex_outer+vy*fc->fy_ex_outer, time);
      fclose (out);
      //printf("disc torque at planet via summation = %lg (x=%lg, y=%lg, ax=%lg, ay=%lg)\n", x*fc->fy_inner-y*fc->fx_inner+x*fc->fy_outer-y*fc->fx_outer,x,y,fc->fx_inner+fc->fx_outer,fc->fy_inner+fc->fy_outer);
    }
  }
}


void UpdateLogSG (psys, outputnb, time)
     PlanetarySystem *psys;
     int outputnb;
     real time;
{
  int k, nb, ip, jp;
  real x, y, rp, tp, vx, vy;
  FILE *out;
  char filename[MAX1D];
  real delta_r, delta_phi;
  extern boolean LogGrid, NGPInterpolation, CICInterpolation, TSCInterpolation;
  real radsgaccatplanet, azisgaccatplanet, sgaccatplanet, xsgaccatplanet, ysgaccatplanet;
  int i, j, myip, myjp, imin, imax, jmin, jmax, myj, nr, ns;
  real myazimuth, dr, dphi, wr, wt;

  nb = psys->nb;
  nr = gr->Nrad;
  ns = gr->Nsec;
  real *radsgacc = gr->Field;
  real *azisgacc = gtheta->Field;
  
  /* Grid spacing along radial and azimuthal directions */
  if (!LogGrid)
    delta_r = Rsup[0]-Rinf[0];        // arithmetic radial spacing
  else
    delta_r = log(Rsup[0]/Rinf[0]);   // logarithmic radial spacing
  delta_phi = AziSup[0]-AziInf[0];
  
  for (k = 0; k < nb; k++) {
    x = psys->x[k];
    y = psys->y[k];
    vx = psys->vx[k];
    vy = psys->vy[k];
    rp = sqrt(x*x+y*y);
    radsgaccatplanet = 0.0;
    azisgaccatplanet = 0.0;

    /* Find out which CPU handles each planet */
    if ( (rp >= Rinf[Zero_or_active]) && (rp < Rsup[Max_or_active-1]) ) {
      /* ----------------------------------------- */
      /* 1) First find in which cell the planet is: (ip,jp) */
      /* ----------------------------------------- */
      ip = Zero_or_active;
      while ( rp >= Rinf[ip] ) ip++;
      ip--;
      tp = atan2(y,x);  // planet's azimuth; beware that atan2 is in [-pi,pi]
      if (tp < AziInf[0])
	tp += (AziSup[NSEC-1]-AziInf[0]);  // add 2pi
      if (tp > AziSup[NSEC-1])
	tp -= (AziSup[NSEC-1]-AziInf[0]);  // subtract 2pi
      if ( (tp < AziInf[0]) || (tp > AziSup[NSEC-1]) )
	printf ("Pb in UpdateLogSG for planet %d: azimuth=%lg\n",k,tp);
      jp = floor(ns*(tp-AziInf[0])/(AziSup[ns-1]-AziInf[0]));
      if (jp == ns)
	jp = 0;
      if (jp == -1) {
	printf ("jp = -1 in interpolation!");
      }

      //printf("x = %lg, y = %lg, rp = %lg , tp = %lg\n", x, y, rp, tp);
      //printf("ip = %d, jp = %d, rinf(ip,ip+1) = %lg %lg, azimuth(jp,jp+1) = %lg %lg\n", ip, jp, Rinf[ip], Rinf[ip+1], AziInf[jp], AziInf[jp+1]);
      
      /* ----------------------------------------- */
      /* 2) Compute interpolated radial and azimuthal SG accelerations
	 at planet's location, knowing that sgacc components are
	 CENTERED IN CELL! */
      /* ----------------------------------------- */
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
	    //printf("dr = %lg, wr = %lg, wt = %lg\n",dr,wr,wt);
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
	    printf ("Pb in step 2 in UpdateLogSG for planet %d, wr=%lg, wt=%lg, i=%d, j=%d\n",k,wr,wt,i,myj);
	  if ( (i >= 0) && (i <= nr-1) ) {
	    /* Now calculate gas SG interpolated accelerations at
	       planet's location */
	    //printf("azisgacc[%d,%d] = %lg, wr=%lg, wt=%lg\n",i,myj,azisgacc[i*ns+myj],wr,wt);
	    radsgaccatplanet += (wr*wt*radsgacc[i*ns+myj]);
	    azisgaccatplanet += (wr*wt*azisgacc[i*ns+myj]);
	  }
	}
      }
      /* --------------------------------------- */
      /* 3) Infer x- and y- SG accelerations at planet location */
      /* --------------------------------------- */
      tp = atan2(y,x);  // planet's azimuth; beware that atan2 is in [-pi,pi]
      xsgaccatplanet = radsgaccatplanet*cos(tp) - azisgaccatplanet*sin(tp);
      ysgaccatplanet = radsgaccatplanet*sin(tp) + azisgaccatplanet*cos(tp);
      
      /* */
      sprintf (filename, "%stqwkSG%d.dat", OUTPUTDIR, k);
      out = fopen (filename, "a");
      if (out == NULL) {
	fprintf (stderr, "Can't open %s\n", filename);
	fprintf (stderr, "Aborted.\n");
      }
      fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb, \
	       x*ysgaccatplanet-y*xsgaccatplanet,			\
	       vx*xsgaccatplanet+vy*ysgaccatplanet,			\
	       radsgaccatplanet, azisgaccatplanet, tp, time);		       
      fclose (out);

      //printf("SG disc torques at planet = %lg (x=%lg, y=%lg, ax=%lg, ay=%lg)\n", x*ysgaccatplanet-y*xsgaccatplanet,x,y,xsgaccatplanet,ysgaccatplanet);
      
    }

  } // end loop over planets
}


void CompareSGandSummationTorques (fc, Rho, psys)
     Force *fc;
     PolarGrid *Rho;
     PlanetarySystem *psys;
{
  int i, j, l, ns;
  real x, y, r, m, smoothing;
  extern boolean SelfGravity;
  
  real *globalforce, *abs, *ord;

  real *torqueSG = torquesg->Field;
  real *azisgacc = gtheta->Field;
  real *torqueSumDisc = torquesumdisc->Field;
  ns = Rho->Nsec;

  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* beware that azisgacc is centred-in-cell so for x and y we
	 need to take the values for cell centres! */
      x = abs[l];
      y = ord[l];
      r = sqrt(x*x+y*y);
      m = 0.0;
      smoothing = compute_smoothing(r);
      ComputeForce (fc, Rho, x, y, smoothing, m, psys);
      globalforce = fc->GlobalForce;
      torqueSumDisc[l] = x*fc->fy_inner-y*fc->fx_inner+x*fc->fy_outer-y*fc->fx_outer;
      if (SelfGravity) {
	torqueSG[l] = r*azisgacc[l];
      }
    }
  }
}
