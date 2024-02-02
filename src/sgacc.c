#include "mp.h"

void compute_sgacc (Rho)
     PolarGrid *Rho;
{
  MPI_Request req1, req3, req4, req5, req6;
  int i, j, l, nr;
  int one_if_odd;
  int stride, mi, mr;
  int ghost_size;
  real normaccr, normacct;
  real *dens, *sgaccr, *sgacctheta;
  dens = Rho->Field;
  nr = Rho->Nrad;
  sgaccr = gr->Field;
  sgacctheta = gtheta->Field;
  stride = 2*(NSEC/2 + 1);
  one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
  
  /* First we compute sg_acc as a convolution product of reduced
     density and kernel arrays. Note that all bufffttabs are
     transposed arrays, since we use flag FFTW_TRANSPOSED_ORDER. */
  SGP_buffft_Accr = (real *) SGP_Accr;
  SGP_buffft_Acct = (real *) SGP_Acct;
  for ( i = 0; i < total_local_size; i++ ) {
    mr = i;   // real part
    mi = i+1; // imaginary part
    SGP_buffft_Accr[mr] = -G*( SGP_buffft_Kr[mr]*SGP_buffft_Sr[mr] -	\
			       SGP_buffft_Kr[mi]*SGP_buffft_Sr[mi] );
    SGP_buffft_Accr[mi] = -G*( SGP_buffft_Kr[mr]*SGP_buffft_Sr[mi] +	\
			       SGP_buffft_Kr[mi]*SGP_buffft_Sr[mr] );
    
    SGP_buffft_Acct[mr] = -G*( SGP_buffft_Kt[mr]*SGP_buffft_St[mr] -	\
			       SGP_buffft_Kt[mi]*SGP_buffft_St[mi] );
    SGP_buffft_Acct[mi] = -G*( SGP_buffft_Kt[mr]*SGP_buffft_St[mi] +	\
			       SGP_buffft_Kt[mi]*SGP_buffft_St[mr] );
    i = i+1;
  }
  rfftwnd_mpi (SGP_fftplan_backward, 1, SGP_Accr, NULL, FFTW_TRANSPOSED_ORDER);
  SGP_buffft_Accr = (real *) SGP_Accr;
  rfftwnd_mpi (SGP_fftplan_backward, 1, SGP_Acct, NULL, FFTW_TRANSPOSED_ORDER);
  SGP_buffft_Acct = (real *) SGP_Acct;
  
  /* The use of argument FFTW_TRANSPOSED_ORDER in above backward
     Fourier Transforms ensures that arrays are not transposed
     anymore. Then, we transfer the exact necessary quantity of
     sg_acceleration arrays from the fftw d.d. to the hydro mesh */
  if ( CPU_Rank != CPU_NoFriend ) {
    if ( CPU_Rank < CPU_Number/2 ) {
      for ( i = 0; i < transfer_size; i++ ) {
	if ( i < transfer_size/2 )
	  ffttohydro_transfer[i] = SGP_buffft_Accr[(ifront+1-CPUOVERLAP)*stride + i];
	else
	  ffttohydro_transfer[i] = SGP_buffft_Acct[(ifront+1-CPUOVERLAP)*stride + i - transfer_size/2];
      }  
      MPI_Isend (ffttohydro_transfer, transfer_size, MPI_DOUBLE, CPU_Friend, 40, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    else {
      MPI_Irecv (ffttohydro_transfer_friend, transfer_size_friend, MPI_DOUBLE, CPU_Friend, 40, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
  }
  
  /* We now compute sg_acceleration arrays on the hydro mesh */
  if ( CPU_Rank < (CPU_Number+one_if_odd)/2 )  {
    if ( CPU_Rank == 0 ) {
      for ( i = 0 ; i < nr; i++ ) {
	for ( j = 0; j < NSEC; j++ ) {
	  l = i*NSEC + j;
	  SG_Accr[l] = SGP_buffft_Accr[i*stride + j];
	  SG_Acct[l] = SGP_buffft_Acct[i*stride + j];
	}
      }
    }
    else {
      for ( i = Zero_or_active ; i < nr; i++ ) {
	for ( j = 0; j < NSEC; j++ ) {
	  l = i*NSEC + j;
	  SG_Accr[l] = SGP_buffft_Accr[(i-Zero_or_active)*stride + j];
	  SG_Acct[l] = SGP_buffft_Acct[(i-Zero_or_active)*stride + j];
	}
      }
    }
  }
  if ( CPU_Rank >= (CPU_Number+one_if_odd)/2 ) {
    if ( CPU_Rank == CPU_Highest ) {
      for ( i = 0 ; i < nr; i++ ) {
	for ( j = 0; j < NSEC; j++ ) {
	  l = i*NSEC + j;
	  SG_Accr[l] = ffttohydro_transfer_friend[i*stride + j];
	  SG_Acct[l] = ffttohydro_transfer_friend[transfer_size_friend/2 + i*stride + j];
	}
      }
    }
    else {
      for ( i = 0 ; i < Max_or_active; i++ ) {
	for ( j = 0; j < NSEC; j++ ) {
	  l = i*NSEC + j;
	  SG_Accr[l] = ffttohydro_transfer_friend[i*stride + j];
	  SG_Acct[l] = ffttohydro_transfer_friend[transfer_size_friend/2 + i*stride + j];
	}
      }
    }
  }
  
  /* Now we exchange the correct amount of sg_acceleration between
     cpus to fill ghosts. */
  ghost_size = CPUOVERLAP * NSEC;
  if ( CPU_Number > 1 ) {
    if ( (CPU_Rank > 0) && (CPU_Rank < (CPU_Number+one_if_odd)/2) ) {
      MPI_Isend (&SG_Acct[Zero_or_active*NSEC], ghost_size, MPI_DOUBLE, CPU_Prev, 60, MPI_COMM_WORLD, &req3);
      MPI_Wait (&req3, &fargostat);
      MPI_Irecv (&SG_Acct[0], ghost_size, MPI_DOUBLE, CPU_Prev, 61, MPI_COMM_WORLD, &req4);
      MPI_Wait (&req4, &fargostat);
    }
    if ( (CPU_Rank >= (CPU_Number+one_if_odd)/2) && (CPU_Rank != CPU_Highest) ) {
      MPI_Irecv (&SG_Acct[Max_or_active*NSEC], ghost_size, MPI_DOUBLE, CPU_Next, 60, MPI_COMM_WORLD, &req3);
      MPI_Wait (&req3, &fargostat);
      MPI_Isend (&SG_Acct[(Max_or_active-CPUOVERLAP)*NSEC], ghost_size, MPI_DOUBLE, CPU_Next, 61, MPI_COMM_WORLD, &req4);
      MPI_Wait (&req4, &fargostat);
    }
    if ( (CPU_Rank > 0) && (CPU_Rank < (CPU_Number+one_if_odd)/2) ) {
      MPI_Isend (&SG_Accr[Zero_or_active*NSEC], ghost_size, MPI_DOUBLE, CPU_Prev, 50, MPI_COMM_WORLD, &req5);
      MPI_Wait (&req5, &fargostat);
      MPI_Irecv (&SG_Accr[0], ghost_size, MPI_DOUBLE, CPU_Prev, 51, MPI_COMM_WORLD, &req6);
      MPI_Wait (&req6, &fargostat);
    }
    if ( (CPU_Rank >= (CPU_Number+one_if_odd)/2) && (CPU_Rank != CPU_Highest) ) {
      MPI_Irecv (&SG_Accr[Max_or_active*NSEC], ghost_size, MPI_DOUBLE, CPU_Next, 50, MPI_COMM_WORLD, &req5);
      MPI_Wait (&req5, &fargostat);
      MPI_Isend (&SG_Accr[(Max_or_active-CPUOVERLAP)*NSEC], ghost_size, MPI_DOUBLE, CPU_Next, 51, MPI_COMM_WORLD, &req6);
      MPI_Wait (&req6, &fargostat);
    }
  }
  
  /* We do not forget to renormalize acc arrays! */
  for ( i = 0 ; i < nr; i++ ) {
    normaccr = SGP_rstep * SGP_tstep / ( (real)(2*GLOBALNRAD) * (real)NSEC );
    normacct = normaccr;
    normaccr /= sqrt(Rmed[i] / GlobalRmed[0]);
    normacct /= ( Rmed[i] / GlobalRmed[0] * sqrt(Rmed[i] / GlobalRmed[0]) );
    for ( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      SG_Accr[l] *= normaccr;
      SG_Acct[l] *= normacct;
    }
  }
  
  /* Eventually, we take the compensation from selfforce into
     account */
  for ( i = 0 ; i < nr; i++ ) {
    for ( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      if ( (i+IMIN) < GLOBALNRAD ) {
	SG_Accr[l] += G*dens[l]*SGP_rstep*SGP_tstep / SGP_eps;
	sgaccr[l] = SG_Accr[l];
	sgacctheta[l] = SG_Acct[l];
      }
    }
  }
}
