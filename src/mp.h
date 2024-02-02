#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "fondam.h"
#include "types.h"
#include "proto.h"
/* PARALLEL LIBS */
#ifdef _PARALLEL
#include <mpi.h>
#ifdef _FFTW
#include <fftw_mpi.h>
#include <rfftw_mpi.h>
#else
#include "fftw_dummy.h"
#endif
#endif
/* SEQUENTIAL LIBS */
#ifndef _PARALLEL
#include "mpi_dummy.h"
#ifdef _FFTW
#include <fftw.h>
#include <rfftw.h>
#else
#include "fftw_dummy.h"
#endif
#endif

#include "sgproto.h"

#ifndef __LOCAL
#include "param.h"
#include "global_ex.h"
#include "sg_ex.h"
#else
#include "param_noex.h"
#include "global.h"
#include "sg.h"
#endif
#include <sys/times.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#ifdef _TRAP_FPE
#include <signal.h>
#include <fenv.h>
#endif
