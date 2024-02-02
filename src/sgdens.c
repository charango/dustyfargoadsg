#include "mp.h"

void compute_fftdensity (Rho)
     PolarGrid *Rho;
{
  MPI_Request req;
  int i, j, l;
  int ih, lh, ir;
  int one_if_odd, stride;
  real *dens;
  dens = Rho->Field;
  stride = 2*(NSEC/2+1);
  
  /* We communicate the hydro density field to the fftw domain
     decomposition (d. d.) */
  one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
  if ( CPU_Rank != CPU_NoFriend ) {
    if ( CPU_Rank >= (CPU_Number+one_if_odd)/2 )
      MPI_Isend (&dens[Zero_or_active*NSEC], active_hydro_totalsize, MPI_DOUBLE, CPU_Friend, 30, MPI_COMM_WORLD, &req);
    else
      MPI_Irecv (&dens_friend[0], active_hydro_totalsize_friend, MPI_DOUBLE, CPU_Friend, 30, MPI_COMM_WORLD, &req);
    MPI_Wait (&req, &fargostat);
  }
  
  if ( (CPU_Rank < CPU_Number/2) && (CPU_Rank!=CPU_NoFriend) ) {
    for ( i = 0; i < ifront+1; i++ ) {
      for ( j = 0; j < NSEC; j++ ) {
	l = i*stride + j;
	ih = i+Zero_or_active;
	lh = ih*NSEC + j;
	SGP_Sr[l] = dens[lh] * sqrt(Rmed[ih] / GlobalRmed[0]);
	SGP_St[l] = SGP_Sr[l] * Rmed[ih] / GlobalRmed[0];
      }
    }
    for ( i = ifront+1; i < local_Nx; i++ ) {
      for ( j = 0; j < NSEC; j++ ) {
	l = i*stride + j;
	ih = i-(ifront+1);
	lh = ih*NSEC + j;
	ir = i+IMIN+Zero_or_active;
	if ( (i+local_i_start) < GLOBALNRAD ) {
	  SGP_Sr[l] = dens_friend[lh] * sqrt(GlobalRmed[ir] / GlobalRmed[0]);
	  SGP_St[l] = SGP_Sr[l] * GlobalRmed[ir] / GlobalRmed[0];
	}
	else {
	  SGP_Sr[l] = 0.;
	  SGP_St[l] = 0.;
	}
      }
    }
  }
  if ( CPU_Rank == CPU_NoFriend ) {
    for ( i = 0; i < local_Nx; i++ ) {
      for ( j = 0; j < NSEC; j++ ) {
	l = i*stride + j;
	ih = i+Zero_or_active;
	if ( (i+local_i_start) < GLOBALNRAD ) {
	  lh = ih*NSEC + j;
	  SGP_Sr[l] = dens[lh] * sqrt(Rmed[ih] / GlobalRmed[0]);
	  SGP_St[l] = SGP_Sr[l] * Rmed[ih] / GlobalRmed[0];
	}
	else {
	  SGP_Sr[l] = 0.;
	  SGP_St[l] = 0.;
	}
      }
    }
  }
  if ( (CPU_Rank >= CPU_Number/2) && (CPU_Rank!=CPU_NoFriend) ) {
    for ( i = 0; i < local_Nx; i++ ) {
      for ( j = 0; j < NSEC; j++ ) {
	l = i*stride + j;
	if ( (i+local_i_start) >= GLOBALNRAD ) {
	  SGP_Sr[l] = 0.;
	  SGP_St[l] = 0.;
	}
      }
    }
  }
  
  /* Now we can compute the in-place ffts of reduced density arrays. */
  rfftwnd_mpi (SGP_fftplan_forward, 1, SGP_Sr, NULL, FFTW_TRANSPOSED_ORDER);
  SGP_buffft_Sr = (real *) SGP_Sr;
  rfftwnd_mpi (SGP_fftplan_forward, 1, SGP_St, NULL, FFTW_TRANSPOSED_ORDER);
  SGP_buffft_St = (real *) SGP_St;
}
