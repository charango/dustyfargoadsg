/** \file Theo.c

A few functions that manipulate the surface density, internal energy
and cooling time profiles.

*/

#include "mp.h"

extern real ScalingFactor; 

/* Surface density */
real Sigma(r)
     real r;
{
  real cavity = 1.0;
  real vr_over_cs, sigmabg;
  real jump, width, F;
  extern boolean TailOffGauss, TailOffIn, TailOffAurelien, TailOffStype, TailOffGI, ExponentialCutoff, CavityTorque;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */
  sigmabg = cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
  if (TailOffGauss) {
    /* We take a Gaussian initial profile */
    sigmabg = SIGMA0*exp(-0.5*pow(r-1.0,2.0)*pow(2.0*ASPECTRATIO,-2.0)) + DENSITYJUMP*SIGMA0;
  }
  if (TailOffIn) {
    /* We take here the prescription in Heemskrek, Papaloizou & Savonije 1992 */
    if (r < 1.4)
      sigmabg *= pow(r-GlobalRmed[0]+0.1,2.0);
  }
  if (ExponentialCutoff) {
    /* We take here the prescription in Owen et al. 2011b */
    //sigmabg *= exp(-r*FACTORUNITLENGTH/CUTDIST); // CUTDIST is in AU in .par file
    sigmabg *= exp(-pow(r*FACTORUNITLENGTH/CUTDIST,EXPONENTIALCUTOFFINDEX)); // CUTDIST is in AU in .par file
  }
  if (TailOffAurelien) {
    /* July 2022 */
    //sigmabg *= exp(-pow(r/4.0,8.0))*exp(-pow(r,-1.2));  // use SIGMA0 = 2.5e-3, SIGMASLOPE = 0.5
    sigmabg *= exp(-pow(r/4.0,8.0))*exp(-pow(r/1.5,-1.1));  // use SIGMA0 = 4.0e-3, SIGMASLOPE = 2.5
  }
  if (TailOffStype) {
    /* July 2023 */
    sigmabg *= exp(-pow(r*FACTORUNITLENGTH/CUTDIST,EXPONENTIALCUTOFFINDEX)); // CUTDIST is in AU in .par file
    sigmabg *= exp(-pow(r/1.0,-1.5));
  }
  if (TailOffGI) {
    /* July 2022 */
    sigmabg *= exp(-pow(r/3.3,8.0))*exp(-pow(r/1.6,-1.1)); // use SIGMA0 = 1.0e-1, SIGMASLOPE = 2.5
  }
  if (CavityTorque) {
    /* June 2022 (Guillaume Robert's internship) */
    vr_over_cs = 0.5*(FRACINT+FRACEXT) + 0.5*(FRACINT-FRACEXT)*tanh((CAVITYRADIUS-r)/CAVITYWIDTH);  // positive definite here!
    sigmabg = SIGMA0 * (FRACEXT/vr_over_cs) * pow(GlobalRmed[NRAD-1]/r,0.5+FLARINGINDEX);
  }
  
  return sigmabg;
}

real DSigma(r)
real r;
{
  real sigmabg;
  extern boolean TailOffGauss, RestartWithNewDust;
  sigmabg = Sigma(r)*DUSTTOGASDENSITYRATIO;
  if (RestartWithNewDust) {
    if ( (r >= RMINDUST) && (r <= RMAXDUST) )
      sigmabg = DSIGMA0*pow(r,-DSIGMASLOPE);
    else {
      if (r < RMINDUST)
	sigmabg =  DSIGMA0*pow(RMINDUST,-DSIGMASLOPE) * exp(-0.5*pow(r-RMINDUST,2.0)*pow(5.0*DAspectRatio(RMINDUST)*pow(RMINDUST,1.0+DFLARINGINDEX),-2.0));
      if (r > RMINDUST)
	sigmabg =  DSIGMA0*pow(RMAXDUST,-DSIGMASLOPE) * exp(-0.5*pow(r-RMAXDUST,2.0)*pow(5.0*DAspectRatio(RMAXDUST)*pow(RMAXDUST,1.0+DFLARINGINDEX),-2.0));
    }
    sigmabg += floordens;
  }
  return sigmabg;
}

void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void FillDSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    DSigmaMed[i] = DSigma(Rmed[i]);
    DSigmaInf[i] = DSigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

void RefillDSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    DSigmaMed[i] = moy;
  }
  DSigmaInf[0] = DSigmaMed[0];
  for (i = 1; i < nr; i++) {
    DSigmaInf[i] = (DSigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   DSigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  real cavity = 1.0;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  }
  else {
    energy0 = R/MU/(ADIABATICINDEX-1.0)*Sigma(r)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);
    //energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
  }
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  return cavity*ScalingFactor*energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Temperature prescription time */
real PrescTime(r)
     real r;
{
  real pt0;
  pt0 = PRESCTIME0*pow(r,2.0+2.0*FLARINGINDEX);
  return pt0;
}

void FillPrescTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    PrescTimeMed[i] = PrescTime(Rmed[i]);
}


/* Beta cooling prescription time */
real ComputeBetaCooling(r)
     real r;
{
  real pt0;
  pt0 = BETACOOLINGTIME*pow(r,-BETACOOLINGSLOPE);
  return pt0;
}
