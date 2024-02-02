/** \file TransportEuler.c

    Functions that handle the transport substep of a hydrodynamical time
    step.  The FARGO algorithm is implemented here. The transport is
    performed in a manner similar to what is done for the ZEUS code (Stone
    & Norman, 1992), except for the momenta transport (we define a left
    and right momentum for each zone, which we declare zone centered; we
    then transport then normally, and deduce the new velocity in each zone
    by a proper averaging).

*/

#include "mp.h"

static int Nshift[MAX1D];
static real *TempShift;
static boolean NoSplitAdvection[MAX1D];
static boolean UniformTransport, DustFluid, ShortFrictionTimeApproximation;

static real *dq, *Ddq;
static PolarGrid *RadMomP, *RadMomM, *ThetaMomP, *ThetaMomM, *ExtLabel, *VthetaRes,  *DRadMomP, *DRadMomM, *DThetaMomP, *DThetaMomM, *DVthetaRes, *QRStarg, *Deri, *Kappa, *KappaStar;
static PolarGrid *Work, *QRStar, *Elongations,  *DWork, *DQRStar;

extern int TimeStep;
extern boolean OpenInner, FastTransport;
extern boolean EnergyEquation;
extern boolean  AdvecteLabel;
real LostMass = 0.0, LostMassd = 0.0;
real AccRate;
real AccRated;

void InitTransport () 
{
  RadMomP      = CreatePolarGrid(NRAD, NSEC, "RadMomP");
  RadMomM      = CreatePolarGrid(NRAD, NSEC, "RadMomM");
  ThetaMomP    = CreatePolarGrid(NRAD, NSEC, "ThetaMomP");
  ThetaMomM    = CreatePolarGrid(NRAD, NSEC, "ThetaMomM");
  Work         = CreatePolarGrid(NRAD, NSEC, "WorkGrid");
  QRStar       = CreatePolarGrid(NRAD, NSEC, "QRStar");
  ExtLabel     = CreatePolarGrid(NRAD, NSEC, "ExtLabel");
  VthetaRes    = CreatePolarGrid(NRAD, NSEC, "VThetaRes");
  Elongations  = CreatePolarGrid(NRAD, NSEC, "Elongations");
  TempShift    = (real *)malloc(NRAD*NSEC*sizeof(real));
  dq           = (real *)malloc(NRAD*NSEC*sizeof(real));
  if (DustFluid) {
    DRadMomP      = CreatePolarGrid(NRAD, NSEC, "DRadMomP");
    DRadMomM      = CreatePolarGrid(NRAD, NSEC, "DRadMomM");
    DThetaMomP    = CreatePolarGrid(NRAD, NSEC, "DThetaMomP");
    DThetaMomM    = CreatePolarGrid(NRAD, NSEC, "DThetaMomM");
    DWork         = CreatePolarGrid(NRAD, NSEC, "DWorkGrid");
    DQRStar       = CreatePolarGrid(NRAD, NSEC, "DQRStar");
    QRStarg       = CreatePolarGrid(NRAD, NSEC, "QRStarg");
    Deri          = CreatePolarGrid(NRAD, NSEC, "Deri");
    Kappa         = CreatePolarGrid(NRAD, NSEC, "Kappa");
    KappaStar     = CreatePolarGrid(NRAD, NSEC, "KappaStar");
    DVthetaRes    = CreatePolarGrid(NRAD, NSEC, "DVThetaRes");
    Ddq           = (real *)malloc(NRAD*NSEC*sizeof(real));
  }
}

/* ============================================ */
/* Function Transport handles the directional
   -splitting transport scheme */
/* ============================================ */
void Transport (Rho, Vrad, Vtheta, Energy, Label, dt)
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
     real dt;
{
  /* First we calculate the four momenta (left/right,
     radial/azimuthal) */
  ComputeLRMomenta (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES)
    ComputeExtQty (Rho, Label, ExtLabel); 
  /* We first transport the four moments, the energy and the density
     in the radial direction */
  OneWindRad (Rho, Vrad, Energy, dt);
  /* We proceed in the same way in the azimuthal direction, except
     that we make use of the FARGO algorithm */
  OneWindTheta (Rho, Vtheta, Energy, dt);
  /* We finally update the velocities with the transported momenta */
  ComputeVelocities (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES) 
    ComputeSpeQty (Rho, Label, ExtLabel);
}

/* ------------------------------- */
/*  Calculation of the LR momenta  */
/* ------------------------------- */
void ComputeLRMomenta (Rho, Vrad, Vtheta)
     PolarGrid *Rho, *Vrad, *Vtheta;
{
  int i,j,l,lip,ljp,nr,ns;
  real *vr,*vt,*rho;
  real *rp, *rm, *tp, *tm;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  rp = RadMomP->Field;
  rm = RadMomM->Field;
  tp = ThetaMomP->Field;
  tm = ThetaMomM->Field; 
  
#pragma omp parallel for private(j,l,lip,ljp)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1)
	ljp = i*ns;
      rp[l] = rho[l]*vr[lip];
      rm[l] = rho[l]*vr[l];
      tp[l] = rho[l]*(vt[ljp]+Rmed[i]*OmegaFrame)*Rmed[i]; /* it is the angular momentum */
      tm[l] = rho[l]*(vt[l]+Rmed[i]*OmegaFrame)*Rmed[i];
    }
  }
}

/* ------------------------------------------ */
/*  Functions regarding the radial transport  */
/* ------------------------------------------ */
void OneWindRad (Rho, Vrad, Energy, dt)
     PolarGrid *Rho, *Vrad, *Energy;
     real dt;
{
  ComputeStarRad (Rho, Vrad, RhoStar, dt);
  ActualiseGas (RhoInt, Rho);
  VanLeerRadial (Vrad, RadMomP, dt); 
  VanLeerRadial (Vrad, RadMomM, dt);  
  VanLeerRadial (Vrad, ThetaMomP, dt);
  VanLeerRadial (Vrad, ThetaMomM, dt);
  if (EnergyEquation)
    VanLeerRadial (Vrad, Energy, dt);
  if (AdvecteLabel == YES)
    VanLeerRadial (Vrad, ExtLabel, dt);
  LostMass += VanLeerRadial (Vrad, Rho, dt); /* MUST be the last line */
}

void ComputeStarRad (Qbase, Vrad, QStar, dt)
     PolarGrid *Qbase, *Vrad, *QStar;
     real dt;
{
  int i,j,l,lq,lip,lim,nr,ns;
  real *qb, *qs, *vr;
  real dqp, dqm;
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qb = Qbase->Field;
  qs = QStar->Field;
  vr = Vrad->Field;
#pragma omp parallel for private (i,l,lq,lip,lim,dqm,dqp)   
  for (j = 0; j < ns; j++) {
    for (i = 0; i < nr; i++) {
      l = j+i*ns;
      lq= i+j*nr;
      lip = l+ns;
      lim = l-ns;
      if ((i == 0) || (i == nr-1)) dq[lq] = 0.0;
      else {
	dqm = (qb[l]-qb[lim])*InvDiffRmed[i];
	dqp = (qb[lip]-qb[l])*InvDiffRmed[i+1];
	if (dqp * dqm > 0.0)
	  dq[lq] = 2.0*dqp*dqm/(dqp+dqm);
	else
	  dq[lq] = 0.0;
      }
    }
    // first index modified from the original public source files
    //for (i = 0; i < nr; i++) {
    for (i = 1; i < nr; i++) {
      l = j+i*ns;
      lq= i+j*nr;
      lip = l+ns;
      lim = l-ns;
      if (vr[l] > 0.0)
	qs[l] = qb[lim]+(Rmed[i]-Rmed[i-1]-vr[l]*dt)*0.5*dq[lq-1];
      else
	qs[l] = qb[l]-(Rmed[i+1]-Rmed[i]+vr[l]*dt)*0.5*dq[lq];
    }
    qs[j] = qs[j+ns*nr] = 0.0;
  }
}

real VanLeerRadial (Vrad, Qbase, dt)
     PolarGrid *Vrad, *Qbase;
     real dt;
{
  int i,j,nr,ns,l,lip;
  real *qrs, *rhos, *vr, *qb;
  real dtheta, varq;
  real LostByDisk=0.0;

  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarRad (Work, Vrad, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vr = Vrad->Field;
  qb=  Qbase->Field;
  dtheta = (PMAX-PMIN)/(real)ns;
  
#pragma omp parallel for private(j,l,lip,varq)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      lip=l+ns;
      varq =dt*dtheta*Rinf[i]*qrs[l]*rhos[l]*vr[l];
      varq-=dt*dtheta*Rsup[i]*qrs[lip]*rhos[lip]*vr[lip];
      qb[l] += varq*InvSurf[i];
      //if ((i == 0) && (OpenInner == YES))
      if (i == 0)
#pragma omp atomic
	LostByDisk += varq;
    }
  }
  return LostByDisk;
}

/* --------------------------------------------- */
/*  Functions regarding the azimuthal transport  */
/* --------------------------------------------- */
void OneWindTheta (Rho, Vtheta, Energy, dt)
     PolarGrid *Rho, *Vtheta, *Energy;
     real dt;
{
  ComputeThetaElongations (Vtheta, dt);
  ComputeAverageThetaVelocities (Vtheta, dt);
  ComputeResiduals (Vtheta, dt);
  ComputeConstantResidual (Vtheta, dt);	
  /* Constant residual is in Vtheta from now on */
  UniformTransport = NO;
  QuantitiesAdvection (Rho, VthetaRes, Energy, dt);
  UniformTransport = YES;
  QuantitiesAdvection (Rho, Vtheta, Energy, dt);
  AdvectSHIFT (RadMomP);
  AdvectSHIFT (RadMomM);
  AdvectSHIFT (ThetaMomP);
  AdvectSHIFT (ThetaMomM);
  if (EnergyEquation)
    AdvectSHIFT (Energy);
  if (AdvecteLabel == YES)
    AdvectSHIFT (ExtLabel);
  AdvectSHIFT (Rho);
}
/* Hereafter are the new specific procedures to the fast algorithm */
void ComputeThetaElongations (Vtheta, dt)
     PolarGrid *Vtheta;
     real dt;
{
  int i,j,l,nr,ns;
  real *vt, *elong;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  vt = Vtheta->Field;
  elong = Elongations->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      elong[l] = vt[l]*dt;
    }
  }
}

void ComputeAverageThetaVelocities (Vtheta, dt)
     PolarGrid *Vtheta;
     real dt;
{
  int i,j,l,nr,ns;
  real *elong;
  real moy, invdt;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  elong = Elongations->Field;
  invdt = 1.0/dt;
#pragma omp parallel for private(j,l,moy)
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      moy += elong[l]*invdt;
    }
    VMed[i] = moy/(real)ns;
  }
}

void ComputeResiduals (Vtheta, dt)
     PolarGrid *Vtheta;
     real dt;
{
  int i,j,l,nr,ns;
  real *vtr, *elong, invdt;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  vtr= VthetaRes->Field;
  elong = Elongations->Field;
  invdt = 1.0/dt;
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      vtr[l] = elong[l]*invdt-VMed[i];
    }
  }
}

void ComputeConstantResidual (Vtheta, dt)
     PolarGrid *Vtheta;
     real dt;
{
  int i,j,l,nr,ns;
  long nitemp;
  real *vt, *vres, Ntilde, Nround, maxfrac, invdt, dpinvns;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  vt = Vtheta->Field;
  vres = VthetaRes->Field;
  invdt = 1.0/dt;
  dpinvns = (PMAX-PMIN)/(real)ns;
  if (FastTransport == YES)
    maxfrac = 1.0;		/* Fast algorithm */
  else
    maxfrac = 0.0;
#pragma omp parallel for private(Ntilde,nitemp,Nround,j,l)
  for (i = 0; i < nr; i++) {
    Ntilde = VMed[i]*InvRmed[i]*dt*(real)ns/(PMAX-PMIN);
    Nround = floor(Ntilde+0.5);
    nitemp = (long)Nround;
    Nshift[i] = (long)nitemp;
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      vt[l] = (Ntilde-Nround)*Rmed[i]*invdt*dpinvns;
    }
    if (maxfrac < 0.5) { 
      NoSplitAdvection[i] = YES;
      for (j = 0; j < ns; j++) {
	l=j+i*ns;
	vres[l] = vt[l]+vres[l];
	vt[l] = 0.0;
      }
    } else {
      NoSplitAdvection[i] = NO;
    }
  }
}

void AdvectSHIFT (array)
     PolarGrid *array;
{
  int i,j,ji,l,li,nr,ns;
  real *val;
  val = array->Field;
  nr  = array->Nrad;
  ns  = array->Nsec;
#pragma omp parallel for private (j,ji,l,li)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      ji = j-Nshift[i];
      while (ji < 0) ji += ns;
      while (ji >= ns) ji -= ns;
      l = j+i*ns;
      li= ji+i*ns;
      TempShift[l]=val[li];
    }
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      val[l] = TempShift[l];
    }
  }
}
/* End of new specific procedures to the fast algorithm */
void QuantitiesAdvection (Rho, Vtheta, Energy, dt)
     PolarGrid *Rho, *Vtheta, *Energy;
     real dt;
{
  ComputeStarTheta (Rho, Vtheta, RhoStar, dt);
  ActualiseGas (RhoInt, Rho);
  VanLeerTheta (Vtheta, RadMomP, dt);
  VanLeerTheta (Vtheta, RadMomM, dt);    
  VanLeerTheta (Vtheta, ThetaMomP, dt);    
  VanLeerTheta (Vtheta, ThetaMomM, dt); 
  if (EnergyEquation)
    VanLeerTheta (Vtheta, Energy, dt); 
  if (AdvecteLabel == YES)
    VanLeerTheta (Vtheta, ExtLabel, dt); 
  VanLeerTheta (Vtheta, Rho, dt); /* MUST be the last line */
}

void ComputeStarTheta (Qbase, Vtheta, QStar, dt)
     PolarGrid *Qbase, *Vtheta, *QStar;
     real dt;
{
  int i,j,l,ljp,ljm,jm,nr,ns;
  real *qb, *qs, *vt;
  real dqp, dqm,dxtheta,ksi,invdxtheta;
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qb = Qbase->Field;
  qs = QStar->Field;
  vt = Vtheta->Field;
#pragma omp parallel for private (dxtheta,invdxtheta,l,j,ljp,ljm,dqm,dqp,jm,ksi)   
  for (i = 0; i < nr; i++) {
    dxtheta = (PMAX-PMIN)/(real)ns*Rmed[i];
    invdxtheta = 1.0/dxtheta;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ljp = l+1;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      if (j == ns-1) ljp = i*ns;
      dqm = (qb[l]-qb[ljm]);
      dqp = (qb[ljp]-qb[l]);
      if (dqp * dqm > 0.0)
	dq[l] = dqp*dqm/(dqp+dqm)*invdxtheta;
      else
	dq[l] = 0.0;
    }
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jm = j-1; 
      if (j == 0) jm = ns-1;
      ljm = jm+i*ns;
      ksi=vt[l]*dt;
      if (ksi > 0.0)
	qs[l] = qb[ljm]+(dxtheta-ksi)*dq[ljm];
      else
	qs[l] = qb[l]-(dxtheta+ksi)*dq[l];
    }
  }
}

void VanLeerTheta (Vtheta, Qbase, dt)
     PolarGrid *Vtheta, *Qbase;
     real dt;
{
  int i,j,nr,ns,l,ljp;
  real *qrs, *rhos, *vt, *qb;
  real dxrad, varq, invsurf;
  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarTheta (Work, Vtheta, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vt = Vtheta->Field;
  qb=  Qbase->Field;
#pragma omp parallel for private(dxrad,invsurf,j,l,ljp,varq)
  for (i = 0; i < nr; i++) {
    dxrad = (Rsup[i]-Rinf[i])*dt;
    invsurf = 1.0/Surf[i];
    if ((UniformTransport == NO) || (NoSplitAdvection[i] == NO)) {
      for (j = 0; j < ns; j++) {
	l=j+i*ns;
	ljp=l+1;
	if (j == ns-1) ljp=i*ns;
	varq  = dxrad*qrs[l]*rhos[l]*vt[l];
	varq -= dxrad*qrs[ljp]*rhos[ljp]*vt[ljp];
	qb[l] += varq*invsurf;
      }
    }
  }
}

/* ----------------------------------------- */
/*  Functions regarding the label advection  */
/* ----------------------------------------- */
void ComputeExtQty (Rho, Label, ExtLabel)
     PolarGrid *Rho, *Label, *ExtLabel;
{
  int i,j,l,nr,ns;
  real *rho, *lab, *extlab;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  lab= Label->Field;
  extlab= ExtLabel->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      extlab[l] = rho[l]*lab[l]; 
      /* compressive flow  if line commentarized
	 extlab[l] = lab[l]; */
    }
  }
}

void ComputeSpeQty (Rho, Label, ExtLabel)
     PolarGrid *Rho, *Label, *ExtLabel;
{
  int i,j,l,nr,ns;
  real *rho, *lab, *extlab;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  lab= Label->Field;
  extlab= ExtLabel->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lab[l] = extlab[l]/rho[l]; 
      /* compressive flow if line commentarized
	 lab[l] = extlab[l]; */
    }
  }
}

/* --------------------------------------- */
/*  Update of the velocitied from momenta  */
/* --------------------------------------- */
void ComputeVelocities (Rho, Vrad, Vtheta)
     PolarGrid *Rho, *Vrad, *Vtheta;
{
  int i,j,l,lim,ljm,nr,ns;
  real *vr,*vt,*rho;
  real *rp, *rm, *tp, *tm;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  rp = RadMomP->Field;
  rm = RadMomM->Field;
  tp = ThetaMomP->Field;
  tm = ThetaMomM->Field;

#pragma omp parallel for private(j,l,lim,ljm)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lim = l-ns;
      ljm = l-1;
      if (j == 0)
	ljm = i*ns+ns-1;
      if (i == 0)
	vr[l] = 0.0;
      else		
	vr[l] = (rp[lim]+rm[l])/(rho[l]+rho[lim]+1e-20);
      vt[l] = (tp[ljm]+tm[l])/(rho[l]+rho[ljm]+1e-20)/Rmed[i]-Rmed[i]*OmegaFrame;
      /* It was the angular momentum */
    }
  }
}


void Transportd (Rho, Rhog, Vrad, Vtheta, Label, dt)
     PolarGrid *Rho, *Rhog, *Vrad, *Vtheta, *Label;
     real dt;
{
  ComputeLRMomenta (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES)
    ComputeExtQty (Rho, Label, ExtLabel);
  /* No-Alternate Directionnal Splitting */
  OneWindRadd (Rho, Rhog, Vrad, dt);
  OneWindThetad (Rho, Rhog, Vtheta, dt);
  ComputeVelocities (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES)
    ComputeSpeQty (Rho, Label, ExtLabel);

}

void OneWindRadd (Rho, Rhog, Vrad, dt) 
     PolarGrid *Rho, *Vrad, *Rhog;
     real dt;
{
  ComputeStarRad (Rho, Vrad, RhoStar, dt);
  ActualiseGas (RhoInt, Rho);
  if(ShortFrictionTimeApproximation == NO){
    VanLeerRadiald (Vrad, RadMomP, dt);
    VanLeerRadiald (Vrad, RadMomM, dt);
    VanLeerRadiald (Vrad, ThetaMomP, dt);
    VanLeerRadiald (Vrad, ThetaMomM, dt);
  }
  if (AdvecteLabel == YES)
    VanLeerRadiald (Vrad, ExtLabel, dt);
  LostMassd += VanLeerRadialdd (Vrad, Rho, Rhog, dt); /* MUST be the last line */
}

void OneWindThetad (Rho, Rhog, Vtheta, dt)
     PolarGrid *Rho, *Vtheta, *Rhog;
     real dt;
{ 
  ComputeThetaElongations (Vtheta, dt);
  ComputeAverageThetaVelocities (Vtheta, dt);
  ComputeResiduals (Vtheta, dt);
  ComputeConstantResidual (Vtheta, dt); /* Constant residual is in Vtheta from now on */
  UniformTransport = NO;
  QuantitiesAdvectiond (Rho, Rhog, VthetaRes, dt);
  UniformTransport = YES;
  QuantitiesAdvectiond (Rho, Rhog, Vtheta, dt);
  AdvectSHIFT (RadMomP);
  AdvectSHIFT (RadMomM);
  AdvectSHIFT (ThetaMomP);
  AdvectSHIFT (ThetaMomM);
  if (AdvecteLabel == YES)
    AdvectSHIFT (ExtLabel);
  AdvectSHIFT (Rho);
}

real VanLeerRadiald (Vrad, Qbase, dt)
     PolarGrid *Vrad, *Qbase;
     real dt;
{ 
  int i,j,nr,ns,l,lip;
  real *qrs, *rhos, *vr, *qb, *work;
  real dtheta, varq;
  real LostByDiskd=0.0;
  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarRad (Work, Vrad, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec; 
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vr = Vrad->Field;
  qb=  Qbase->Field;
  work=Work->Field;
  l=ns;
  dtheta = 2.0*PI/(real)ns;
#pragma omp parallel for private(j,l,lip,varq)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      lip=l+ns;
      varq =dt*dtheta*Rinf[i]*qrs[l]*rhos[l]*vr[l];
      varq-=dt*dtheta*Rsup[i]*qrs[lip]*rhos[lip]*vr[lip];
      qb[l] += varq*InvSurf[i];
    }
  }
  return LostByDiskd;
} 


real VanLeerRadialdd (Vrad, Qbase, Qbaseg, dt)
     PolarGrid *Vrad, *Qbase, *Qbaseg;
     real dt;
{
  int i,j,nr,ns,l,lip,lim;
  real *qrs, *rhos, *vr, *qb, *qbg, *qrsg, *der, *kap, *kaps;
  real dtheta, varq, vis;
  real LostByDiskd=0.0;
  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarRad (Work, Vrad, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vr = Vrad->Field;
  qb=  Qbase->Field;

  qbg= Qbaseg->Field;
  qrsg=QRStarg->Field;
  der= Deri->Field;
  kap=Kappa->Field;
  kaps=KappaStar->Field;
  
  dtheta = 2.0*PI/(real)ns;
  AccRated=0.;
#pragma omp parallel for private(j,l,lip,varq)
  /*
    for (i = 0; i < nr; i++) {
    vis=FViscosity (Rmed[i]); 
    for (j = 0; j < ns; j++) {
    l=j+i*ns;
    lip=l+ns;
    lim=l-ns;
    if (i>0) der[l]=-(qb[l]/qbg[l]-qb[lim]/qbg[lim])* InvDiffRmed[i];
    else der[l]=0.;
    kap[l]=qbg[l]*vis;
    }
    }
    ComputeStarRad (Qbaseg, Deri, QRStarg, dt);
    ComputeStarRad (Kappa, Deri, KappaStar, dt);
  */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      lip=l+ns;
      varq =dt*dtheta*Rinf[i]*(qrs[l]*rhos[l]*vr[l]);
      /*
	if(DustDiff == YES&&Rmed[i]>DUSTDIFFINN){ 
        varq += dt*dtheta*Rinf[i]*kaps[l]*der[l]*qrsg[l];
	}
      */
      varq-=dt*dtheta*Rsup[i]*(qrs[lip]*rhos[lip]*vr[lip]);
      /*
	if(DustDiff == YES&&Rmed[i]>DUSTDIFFINN){
	varq -= dt*dtheta*Rsup[i]*kaps[lip]*der[lip]*qrsg[lip];
	}
      */
      qb[l] += varq*InvSurf[i];
      if ((i == 0) && (OpenInner == YES))
#pragma omp atomic
        LostByDiskd += varq;
      if ((i == 0) && (OpenInner == YES)&& (CPU_Rank==0))
	{       AccRated += dtheta*Rsup[i]*qrs[lip]*rhos[lip]*vr[lip];
	}
    }
  }
  return LostByDiskd;
}

void VanLeerThetad (Vtheta, Qbase, Qbaseg, dt)
     PolarGrid *Vtheta, *Qbase, *Qbaseg;
     real dt;
{
  int i,j,nr,ns,l,ljp, ljm;
  real *qrs, *rhos, *vt, *qb, dxtheta, invdxtheta,*qbg,*qrsg,*der,*kap,*kaps;
  real dxrad, varq, invsurf, vis, dtheta;
  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarTheta (Work, Vtheta, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vt = Vtheta->Field;
  qb=  Qbase->Field;

  qbg= Qbaseg->Field;
  qrsg=QRStarg->Field;
  der= Deri->Field;
  kap=Kappa->Field;
  kaps=KappaStar->Field;

  dtheta = 2.0*PI/(real)ns;
#pragma omp parallel for private(j,l,lip,varq)
  /*
    for (i = 0; i < nr; i++) {
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    invdxtheta = 1.0/dxtheta;
    vis=FViscosity (Rmed[i]); 
    for (j = 0; j < ns; j++) {
    l=j+i*ns;
    ljp=l+1;
    ljm=l-1;
    if (j>0) der[l]=-(qb[l]/qbg[l]-qb[ljm]/qbg[ljm])*invdxtheta;
    else{
    ljm=(ns-1)+i*ns;
    der[l]=-(qb[l]/qbg[l]-qb[ljm]/qbg[ljm])*invdxtheta;
    }
    kap[l]=qbg[l]*vis;
    }
    }
    ComputeStarTheta (Qbaseg, Deri, QRStarg, dt);
    ComputeStarTheta (Kappa, Deri, KappaStar, dt);
  */

#pragma omp parallel for private(dxrad,invsurf,j,l,ljp,varq)
  for (i = 0; i < nr; i++) {
    dxrad = (Rsup[i]-Rinf[i])*dt;
    invsurf = 1.0/Surf[i];
    if ((UniformTransport == NO) || (NoSplitAdvection[i] == NO)) {
      for (j = 0; j < ns; j++) {
        l=j+i*ns;
        ljp=l+1;
        if (j == ns-1) ljp=i*ns;
        varq  = dxrad*(qrs[l]*rhos[l]*vt[l]);
	/*
	  if(DustDiff == YES&&Rmed[i]>DUSTDIFFINN){
	  varq += dxrad*kaps[l]*qrsg[l]*der[l];
	  }
	*/
        varq -= dxrad*(qrs[ljp]*rhos[ljp]*vt[ljp]);
	/*
	  if(DustDiff == YES&&Rmed[i]>DUSTDIFFINN){
	  varq -= dxrad*kaps[ljp]*qrsg[ljp]*der[ljp];
	  }
	*/
        qb[l] += varq*invsurf;
      }
    }
  }
}

void QuantitiesAdvectiond (Rho, Rhog, Vtheta, dt)
PolarGrid *Rho, *Vtheta, *Rhog;
real dt;
{
  ComputeStarTheta (Rho, Vtheta, RhoStar, dt);
  ActualiseGas (RhoInt, Rho);
  if(ShortFrictionTimeApproximation == NO){
    VanLeerTheta (Vtheta, RadMomP, dt);
    VanLeerTheta (Vtheta, RadMomM, dt);
    VanLeerTheta (Vtheta, ThetaMomP, dt);
    VanLeerTheta (Vtheta, ThetaMomM, dt);
  }
  if (AdvecteLabel == YES)
    VanLeerTheta (Vtheta, ExtLabel, dt);
  VanLeerThetad (Vtheta, Rho, Rhog, dt); /* MUST be the last line */
}
