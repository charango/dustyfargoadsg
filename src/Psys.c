/** \file Psys.c

Contains the functions that set up the planetary system configuration.
In addition, the last two functions allow to track the first planet
(number 0) of the planetary system, in order to perform a calculation
in the frame corotating either with this planet or with its
guiding-center.

*/

#include "mp.h"

static real Xplanet, Yplanet;
extern boolean GuidingCenter, BinaryCenter, RetrogradeBinary, ForcedCircular;

int FindNumberOfPlanets (filename)
     char *filename;
{
  FILE *input;
  char s[512];
  int Counter=0;
  input = fopen (filename, "r");
  if (input == NULL) {
    fprintf (stderr, "Error : can't find '%s'.\n", filename);
    prs_exit (1);
  }
  while (fgets(s, 510, input) != NULL) {
    if (isalpha(s[0]))
      Counter++;
  }
  fclose (input);
  return Counter;
}

PlanetarySystem *AllocPlanetSystem (nb)
int nb;
{
  real *mass, *x, *y, *vx, *vy, *acc;
  boolean *feeldisk, *feelothers, *binary;
  int i;
  PlanetarySystem *sys;
  sys  = (PlanetarySystem *)malloc (sizeof(PlanetarySystem));
  if (sys == NULL) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  x    = (real *)malloc (sizeof(real)*(nb+1));
  y    = (real *)malloc (sizeof(real)*(nb+1));
  vy   = (real *)malloc (sizeof(real)*(nb+1));
  vx   = (real *)malloc (sizeof(real)*(nb+1));
  mass = (real *)malloc (sizeof(real)*(nb+1));
  acc  = (real *)malloc (sizeof(real)*(nb+1));
  if ((x == NULL) || (y == NULL) || (vx == NULL) || (vy == NULL) || (acc == NULL) || (mass == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  feeldisk   = (boolean *)malloc (sizeof(real)*(nb+1));
  feelothers = (boolean *)malloc (sizeof(real)*(nb+1));
  binary     = (boolean *)malloc (sizeof(real)*(nb+1));
  if ((feeldisk == NULL) || (feelothers == NULL) || (binary == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  sys->x = x;
  sys->y = y;
  sys->vx= vx;
  sys->vy= vy;
  sys->acc=acc;
  sys->mass = mass;
  sys->FeelDisk = feeldisk;
  sys->FeelOthers = feelothers;
  sys->Binary = binary;
  for (i = 0; i < nb; i++) {
    x[i] = y[i] = vx[i] = vy[i] = mass[i] = acc[i] = 0.0;
    feeldisk[i] = feelothers[i] = YES;
    binary[i] = NO;
  }
  return sys;
}

void FreePlanetary (sys)
     PlanetarySystem *sys;
{
  free (sys->x);
  free (sys->vx);
  free (sys->y);
  free (sys->vy);
  free (sys->mass);
  free (sys->acc);
  free (sys->FeelOthers);
  free (sys->FeelDisk);
  free (sys->Binary);
  free (sys);
}

PlanetarySystem *InitPlanetarySystem (filename)
char *filename;
{
  extern boolean CICPlanet;
  FILE *input;
  char s[512], nm[512], test1[512], test2[512], test3[512], *s1;
  PlanetarySystem *sys;
  int i=0, j, nb;
  float mass, dist, accret;
  real m0, m1, mbin;
  real massinvelocity;
  boolean feeldis, feelothers, binary;
  nb = FindNumberOfPlanets (filename);
  if (CPU_Master)
    printf ("%d planet(s) found.\n", nb);
  sys = AllocPlanetSystem (nb);
  input = fopen (filename, "r");
  sys->nb = nb;
  while (fgets(s, 510, input) != NULL) {
    sscanf(s, "%s ", nm);
    if (isalpha(s[0])) {
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %f %s %s %s", &dist, &mass, &accret, test1, test2, test3);
      if ( CICPlanet ) {
	/* initialization puts planet at the interface between two
	   cells (with eccentricity = 0 only) */
	j = 0;
	while ( GlobalRmed[j] < dist ) j++;
	dist = Radii[j+1];
      }
      sys->mass[i] = (real)mass;
      InitialPlanetMass[i] = 0.0;  // except a restart is done
      FinalPlanetMass[i] = sys->mass[i]; /* Mass of planet i at the
					    end of the calculation */
      // CB: March 2022 only apply MassTaper in calculation of disc
      // potential in PframeForce.c, like in Fargo3D
      if (MASSTAPER > 1e-3)
	      massinvelocity = 0.0;
      else
	      massinvelocity = FinalPlanetMass[i];
      //massinvelocity = FinalPlanetMass[i];
      /* This is done only when MassTaper does not vanish for the
	 first output of the planet's mass in the planet[i].dat
	 files */
      feeldis = feelothers = YES;
      binary = NO;
      if (tolower(*test1) == 'n') feeldis = NO;
      if (tolower(*test2) == 'n') feelothers = NO;
      /* NEW July 2014 */
      if ( (feeldis == YES) && (ForcedCircular) ) {
	mastererr("You can't set FeelDisk to Yes in the .cfg file and ForcedCircular to Yes in the .par file. Aborted!\n");
	prs_exit(1);
      }
      if ( (feelothers == YES) && (ForcedCircular) ) {
	mastererr("You can't set FeelOthers to Yes in the .cfg file and ForcedCircular to Yes in the .par file. Aborted!\n");
	prs_exit(1);
      }
      if (tolower(*test3) == 'y') binary = YES;
      sys->x[i] = (real)dist*(1.0+ECCENTRICITY);
      sys->y[i] = 0.0;
      /*
      sys->vy[i] = (real)sqrt(G*(1.0+sys->mass[i])/dist)*	\
	sqrt( (1.0-ECCENTRICITY)/(1.0+ECCENTRICITY) );
      */
      sys->vy[i] = (real)sqrt(G*(1.0+massinvelocity)/dist)*	\
	sqrt( (1.0-ECCENTRICITY)/(1.0+ECCENTRICITY) );
      sys->vx[i] = -0.0000000001*sys->vy[i];
      sys->acc[i] = accret;
      sys->FeelDisk[i] = feeldis;
      sys->FeelOthers[i] = feelothers;
      sys->Binary[i] = binary;
      i++;
    }
  }
  for (i = 0; i < nb; i++) {
    if (sys->Binary[i] == YES) {
      /* Binary initialization: planet feels binary planet  */
      m0 = sys->mass[0];
      m1 = sys->mass[1];
      mbin = m0+m1;
      if ( i==0 ) {
	sys->x[i] = (real)dist;
	sys->y[i] = (m1/mbin)*BINARYSEPARATION*(1.0+BINARYECCENTRICITY);
	sys->vy[i] = (real)sqrt(G*(1.0+mbin)/sys->x[i]);
	sys->vx[i] = -(real)sqrt(G*m1*m1/mbin/BINARYSEPARATION*		\
				 (1.0-BINARYECCENTRICITY)/(1.0+BINARYECCENTRICITY));
	if (RetrogradeBinary) sys->vx[i] = -sys->vx[i];
      }
      if ( i==1 ) {
	sys->x[i] = sys->x[0];
	sys->y[i] = -(m0/mbin)*BINARYSEPARATION*(1.0+BINARYECCENTRICITY);
	sys->vy[i] = sys->vy[0];
	sys->vx[i] = (real)sqrt(G*m0*m0/mbin/BINARYSEPARATION*		\
				(1.0-BINARYECCENTRICITY)/(1.0+BINARYECCENTRICITY));
	if (RetrogradeBinary) sys->vx[i] = -sys->vx[i];
      }
    }
  }
  HillRadius = sys->x[0] * pow( sys->mass[0]/3., 1./3. );
  return sys;
}

void ListPlanets (sys)
     PlanetarySystem *sys;
{
  int nb;
  int i;
  nb = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < nb; i++) {
    printf ("Planet number %d\n", i);
    printf ("---------------\n");
    printf ("x = %.10f\ty = %.10f\n", sys->x[i],sys->y[i]);
    printf ("vx = %.10f\tvy = %.10f\n", sys->vx[i],sys->vy[i]);
    if (sys->acc[i] == 0.0)
      printf ("Non-accreting.\n");
    else
      printf ("accretion time = %.10f\n", 1.0/(sys->acc[i]));
    if (sys->FeelDisk[i] == YES) {
      printf ("Feels the disk potential\n");
    } else {
      printf ("Doesn't feel the disk potential\n");
    }
    if (sys->FeelOthers[i] == YES) {
      printf ("Feels the other planets potential\n");
    } else {
      printf ("Doesn't feel the other planets potential\n");
    }
    if (sys->Binary[i] == YES) {
      printf ("Is in a binary system\n");
    } else {
      printf ("Is not in a binary system\n");
    }
    printf ("\n");
  }
}

real GetPsysInfo (sys, action)
     PlanetarySystem *sys;
     boolean action;
{
  extern boolean SelfGravity, SGZeroMode;
  extern boolean CorotateWithOuterPlanet;
  real d1, d2, cross;
  real x,y, vx, vy, m, h, d, Ax, Ay, e, a, E, M;
  real xc, yc, vxc, vyc, omega;
  real x0, x1, y0, y1, vx0, vx1, vy0, vy1, m0, m1;
  real arg, PerihelionPA;
  real ri, rip1, dr, sgacc;
  int ipl, nb;
  if (!CorotateWithOuterPlanet) {
    xc = x = sys->x[0];
    yc = y = sys->y[0];
    vxc = vx= sys->vx[0];
    vyc = vy= sys->vy[0];
    m = sys->mass[0]+1.;
  } else {
    nb = sys->nb;
    xc = x = sys->x[nb-1];
    yc = y = sys->y[nb-1];
    vxc = vx= sys->vx[nb-1];
    vyc = vy= sys->vy[nb-1];
    m = sys->mass[nb-1]+1.;
  }
  if (BinaryCenter) {
    x0 = sys->x[0];
    x1 = sys->x[1];
    y0 = sys->y[0];
    y1 = sys->y[1];
    vx0 = sys->vx[0];
    vx1 = sys->vx[1];
    vy0 = sys->vy[0];
    vy1 = sys->vy[1];
    m0 = sys->mass[0];
    m1 = sys->mass[1];
    if (m0 == 0.0)
      m0 = FinalPlanetMass[0];
    if (m1 == 0.0)
      m1 = FinalPlanetMass[1];
    xc = x = (m0*x0 + m1*x1) / (m0+m1);
    yc = y = (m0*y0 + m1*y1) / (m0+m1);
    vxc = vx = (m0*vx0 + m1*vx1) / (m0+m1);
    vyc = vy = (m0*vy0 + m1*vy1) / (m0+m1);
    m = m0+m1+1.;
  }
  h = x*vy-y*vx;
  d = sqrt(x*x+y*y);
  Ax = x*vy*vy-y*vx*vy-G*m*x/d;
  Ay = y*vx*vx-x*vx*vy-G*m*y/d;
  e = sqrt(Ax*Ax+Ay*Ay)/m;
  a = h*h/G/m/(1.-e*e);
  if (e == 0.0) {
    arg = 1.0;
  } else {
    arg = (1.0-d/a)/e;
  }
  if (fabs(arg) >= 1.0) 
    E = PI*(1.-arg/fabs(arg))/2.;
  else
    E = acos((1.0-d/a)/e);
  if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) E= -E;
  M = E-e*sin(E);
  omega = sqrt(m/a/a/a);
  /* Here omega is modified to include self-gravity */
  if ( SelfGravity ) {
    if ( !SGZeroMode )
      mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
    else
      GLOBAL_AxiSGAccr = SG_Accr;
    ipl = 0;
    while (GlobalRmed[ipl] <= a) ipl++;
    ri = GlobalRmed[ipl];
    rip1 = GlobalRmed[ipl+1];
    dr = rip1 - ri;
    sgacc = (a - ri)*GLOBAL_AxiSGAccr[ipl+1] + (rip1 - a)*GLOBAL_AxiSGAccr[ipl];
    sgacc /= dr;
    omega *= (real)sqrt(1. - a*a*sgacc/m);
  }
  PerihelionPA=atan2(Ay,Ax);
  if (GuidingCenter == YES) {
    xc = a*cos(M+PerihelionPA);
    yc = a*sin(M+PerihelionPA);
    vxc = -a*omega*sin(M+PerihelionPA);
    vyc =  a*omega*cos(M+PerihelionPA);
  } 
  if (e < 1e-8) {
    xc = x;
    yc = y;
    vxc = vx;
    vyc = vy;
  }
  switch (action) {
  case MARK: 
    Xplanet = xc;
    Yplanet = yc;
    return 0.;
    break;
  case GET:
    x = xc;
    y = yc;
    vx = vxc;
    vy = vyc;
    d2 = sqrt(x*x+y*y);
    d1 = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
    cross = Xplanet*y-x*Yplanet;
    Xplanet = x;
    Yplanet = y;
    return asin(cross/(d1*d2));
    break;
  case FREQUENCY:
    return omega;
    break;
  }
  return 0.0;
}

void RotatePsys (sys, angle)	/* Rotate by angle '-angle' */
     PlanetarySystem *sys;
     real angle;
{
  int nb;
  int i;
  real sint, cost, xt, yt;
  nb = sys->nb;
  sint = sin(angle);
  cost = cos(angle);
  for (i = 0; i < nb; i++) {
    xt = sys->x[i];
    yt = sys->y[i];
    sys->x[i] = xt*cost+yt*sint;
    sys->y[i] = -xt*sint+yt*cost;
    xt = sys->vx[i];
    yt = sys->vy[i];
    sys->vx[i] = xt*cost+yt*sint;
    sys->vy[i] = -xt*sint+yt*cost;
  }
}
