/** \file types.h

Definition of the structures used in the FARGO code.

*/

#include <sys/times.h>

typedef int     boolean;
typedef double	real;

struct torque_set {
  real          InnerDisk;
  real          OuterDisk;
  real		ExcludeHill;
  real          Total;
};

typedef struct torque_set TorqueSet;

struct triplet {
  real            x;
  real            y;
  real		z;
};

typedef struct triplet Triplet;

struct pair {
  real            x;
  real            y;
};

typedef struct pair Pair;

struct force {
  real fx_inner;
  real fy_inner;
  real fx_ex_inner;
  real fy_ex_inner;
  real fx_outer;
  real fy_outer;
  real fx_ex_outer;
  real fy_ex_outer;
  real *GlobalForce;
};

typedef struct force Force;

struct polargrid {
  int             Nrad;
  int             Nsec;
  real           *Field;
  char           *Name;
};

typedef struct polargrid PolarGrid;

#define		YES	1
#define		NO	0
#define		REAL	1
#define		INT	0
#define		STRING  2
#define 	SINE	0
#define		COSINE	1
#define		ABSCISSA	0
#define		ORDINATE	1
#define		HEIGHT		2
#define		INF		0
#define 	SUP		1
#define         GET             0
#define         MARK            1
#define         FREQUENCY       2
#define         COM_DENSITY     0
#define         COM_VRAD        1
#define         COM_VTHETA      2

#define		MAX1D	16384

struct param {
  char name[80];
  int  type;
  char *variable;
  int read;
  int necessary;
};

typedef struct param Param;

struct timeprocess {
  char name[80];
  clock_t clicks;
};

typedef struct timeprocess TimeProcess;

struct planetary_system {
  int nb;			/* Number of planets */
  real *mass;			/* their masses */
  real *x;			/* their coordinates */
  real *y;
  real *vx;			/* their velocities */
  real *vy;
  real *acc;			/* Their accretion times^-1 */
  char **name;			/* their names */
  boolean *FeelDisk;		/* will "migrate" ? */
  boolean *FeelOthers;		/* will feel other planets ? */
  boolean *Binary;		/* are planets in a binary system ? */
};

typedef struct planetary_system PlanetarySystem;
 

struct DustParticlesSystem {
  real *dustsize;                /* particles size */
  real *dustsize_init;           /* particles size at beignning of simulation */
  real *r;                       /* their radial position */
  real *th;                      /* their azimuthal position */
  real *vr;                      /* their radial velocity */
  real *vth;                     /* their azimuthal velocity */
  real *l;                       /* their specific angular momentum l = rxvth */
  real *dvr;                     /* vr - vr_gas at particle's location */
  real *dvth;                    /* vth - vth_gas at particle's location */
  real *rad_geff_acc;            /* the radial component of their effective gravity acceleration */
  real *azi_geff_acc;            /* the azimuthal component of their effective gravity acceleration */
  real *rad_gradp;               /* the radial component of the gas pressure gradient at particle's location */
  real *azi_gradp;               /* the azimuthal component of the gas pressure gradient at particle's location */
  real *rad_drag_acc;            /* the radial component of their gas drag acceleration */
  real *azi_drag_acc;            /* the azimuthal component of their gas drag acceleration */
  real *rad_ind_acc;             /* the radial component of their indirect acceleration */
  real *azi_ind_acc;             /* the azimuthal component of their indirect acceleration */
  real *rad_sg_acc;              /* the radial component of the gas self-gravitating acceleration */
  real *azi_sg_acc;              /* the azimuthal component of the gas self-gravitating acceleration */
  real *jacobi;                  /* their Jacobi constant */
  real *stokesnb;                /* their stokes number */
  real *gasvratpc;               /* radial gas velocity at pc position */
  real *gasvtatpc;               /* azimuthal gas velocity at pc position */
  real *gasdensatpc;             /* gas density at pc position */
  real *gascsatpc;               /* gas sound speed at pc position */
  real *senddusttoprevCPU;       /* array with information on dust passed to previous CPU */
  real *senddusttonextCPU;       /* array with information on dust passed to next CPU */
  real *recvdustfromprevCPU;     /* array with information on dust obtained from previous CPU */
  real *recvdustfromnextCPU;     /* array with information on dust obtained from next CPU */
};

typedef struct DustParticlesSystem DustSystem;


struct Jacobi {
  real *time;
  real *J;
};

typedef struct Jacobi Jacobi;
