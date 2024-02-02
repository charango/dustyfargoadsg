int CPU_Rank;
int CPU_Number;
boolean CPU_Master;
int CPU_Next, CPU_Prev, CPU_Highest;
real unit_mass, unit_length, unit_temperature, unit_time, mmw, sigma_SB;
int DimToNext, DimToPrev, DimFromPrev, DimFromNext, intfoo = 0;
int dimfxy = 2;
/* ------------------------------------- */
/* Variables specific to fftw mesh split */
/* ------------------------------------- */
int CPU_Friend, CPU_NoFriend;
real *dens_friend;
real *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
real *ffttohydro_transfer, *ffttohydro_transfer_friend;
int local_Nx, local_i_start, local_i_start_friend, total_local_size_friend, local_Nx_friend;
int local_Ny_after_transpose, local_j_start_after_transpose;
int total_local_size, ifront, Zero_or_active_friend;
int transfer_size, transfer_size_friend;
int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;  
/* ------------------------------------- */
int IMIN;
int IMAX;
int Zero_or_active;
int Max_or_active;
int One_or_active;
int MaxMO_or_active;		/* MO: Minus One */
int GLOBALNRAD;
int ievaporation = 0;
real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D], Azimuth[MAX1D], CosAzimuth[MAX1D], SinAzimuth[MAX1D], AziInf[MAX1D], AziSup[MAX1D];
real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper, DustMassTaper;
real DSigmaMed[MAX1D], DSigmaInf[MAX1D];
real EnergyMed[MAX1D], PrescTimeMed[MAX1D];
real InitialPlanetMass[MAX1D], FinalPlanetMass[MAX1D];
real VMed[MAX1D];
real GLOBAL_SoundSpeed[MAX1D], GLOBAL_DustSoundSpeed[MAX1D];
real OmegaFrame, PhysicalTime=0.0, PhysicalTimeInitial;
int TimeStep=0;
real HillRadius, mdcp, mdcp0, exces_mdcp;
real GLOBAL_bufarray[MAX1D];
real Particles_Mass, Particles_Mass_Initial;
boolean Merge, FakeSequential, MonitorIntegral, debug, OnlyInit;
boolean	GotoNextOutput, StoreSigma, StoreEnergy, ViscosityAlpha, DViscosityAlpha, RocheSmoothing;
boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
MPI_Status fargostat;
PolarGrid *CellAbscissa, *CellOrdinate;
PolarGrid *RhoStar, *RhoInt, *Potential, *IndPotential, *TurbPotential, *Pressure, *SoundSpeed, *Temperature, *RadIndAcc, *AziIndAcc, *RadSGAcc, *AziSGAcc, *RadFBAcc, *AziFBAcc, *FBedot, *RadGradP, *AziGradP, *RadiativeKCoeff, *RadiativeChiCoeff;
PolarGrid *DRhoStar, *DRhoInt, *Stokes, *Diag1, *Diag2, *Diag3, *Fdiffrp, *Fdifftp, *DPressure, *DSoundSpeed;
PolarGrid *DivergenceVelocity, *TAURR, *TAUPP, *TAURP, *ViscHeat, *EntropyDiff, *RadiativeDiff, *ThermCool, *Opacity;
PolarGrid *DDivergenceVelocity, *DTAURR, *DTAUPP, *DTAURP;
PolarGrid *gr, *gtheta, *torquesg, *torquesumdisc;
PolarGrid *Test;
PolarGrid *Residual_Rho;
PolarGrid *a_SORarray, *b_SORarray, *c_SORarray, *d_SORarray, *e_SORarray, *f_SORarray;
PolarGrid *Global_a_SORarray, *Global_b_SORarray, *Global_c_SORarray, *Global_d_SORarray, *Global_e_SORarray, *Global_f_SORarray, *Global_tempint;
real *Radii1D, *Rmed1D, *Rsup1D, *Rinf1D, *Sigma1D, *Vrad1D, *cs1D, *viscosity1D, *term1, *term2, *term3;
boolean LogGrid;
boolean OverridesOutputdir;
char NewOutputdir[1024];
real *potturb;
real *GLOBAL_Axidens_Evap;
real VthetaMed[MAX1D], VradMed[MAX1D];
real DVthetaMed[MAX1D], DVradMed[MAX1D];
real Minimum_Stopping_Time;

