/** \file split.c

Split (radially) the mesh among the different processors. When the
disc self-gravity is not included, a simple Round Robin algorithm is
used to achieve a proper load balancing, assuming the cluster to be
homogeneous.  A ring described by a given process has its 5
(CPUOVERLAP) innermost zones that are ghost zones which are filled by
communications with the previous inner process (unless it is the
innermost process itself), and its 5 (CPUOVERLAP) outermost zones that
are ghost zones which are filled by communications with the next outer
process (unless it is the outermost process itself). The "active" part
of the submesh described by a given process is the part of this mesh
that excludes the ghost zones.  For each process, the (local) radial
index of the first active ring is Zero_or_active, and the (local)
radial index of the last active ring is Max_or_active. MaxMO_or_active
is Max_or_active for all processes, except the very last one (the
outermost) for which it is Max_or_active-1 (MO stands for 'Minus
One'). When the disc self-gravity is accounted for, the domain
decomposition is imposed by the fftw library.

*/

#include "mp.h"

void SplitDomain () {
  extern boolean SelfGravity, SGZeroMode;
  /* Variables specific to standard mesh split */
  int remainder;
  int size_low, size_high;
  /* Variables specific to fftw mesh split */
  int IMAX_friend;
  int one_if_odd;
  int sg_split[5], hydro_split[5];
  
  GLOBALNRAD = NRAD;
  /* ----------------------------- */
  /* Standard domain decomposition */
  /* ----------------------------- */
  if ( !SelfGravity || SGZeroMode ) {
    size_low = NRAD/CPU_Number;
    size_high = size_low+1;
    remainder = NRAD % CPU_Number;
    if (size_low < 2*CPUOVERLAP) {
      mastererr("The number of processes is too large\n");
      mastererr("or the mesh is radially too narrow.\n");
      prs_exit(1);
    }
    if (CPU_Rank < remainder) {
      IMIN = size_high*CPU_Rank;
      IMAX = IMIN+size_high-1;
    } else {
      IMIN = size_high*remainder+(CPU_Rank-remainder)*size_low;
      IMAX = IMIN+size_low-1;
    }
    if (CPU_Rank > 0) IMIN -= CPUOVERLAP;
    if (CPU_Rank < CPU_Number-1) IMAX +=CPUOVERLAP;
    NRAD = IMAX-IMIN+1;
    Zero_or_active = CPUOVERLAP * (CPU_Rank > 0 ? 1 : 0);
    One_or_active = 1+(CPUOVERLAP-1) * (CPU_Rank > 0 ? 1 : 0);
    Max_or_active = NRAD-CPUOVERLAP*(CPU_Rank < CPU_Number-1 ? 1 : 0);
    MaxMO_or_active = NRAD-1-(CPUOVERLAP-1)*(CPU_Rank < CPU_Number-1 ? 1 : 0);
    CPU_Highest = CPU_Number-1;
    if ( CPU_Rank > 0 )
      CPU_Prev = CPU_Rank-1;
    if ( CPU_Rank < CPU_Highest )
      CPU_Next = CPU_Rank+1;
  }
  /* ---------------------------------------- */
  /* Domain decomposition imposed by fftw_mpi */
  /* ---------------------------------------- */
  if ( SelfGravity && !SGZeroMode ) {
    /* Here we create plan to calculate forward ffts */
    SGP_fftplan_forward = rfftw2d_mpi_create_plan(MPI_COMM_WORLD,
						  2*GLOBALNRAD, NSEC, 
						  FFTW_REAL_TO_COMPLEX, 
						  FFTW_MEASURE);
    /* It is this function that constrains the domain decomposition in
       the fftw mesh */
    rfftwnd_mpi_local_sizes(SGP_fftplan_forward, &local_Nx, &local_i_start,
			    &local_Ny_after_transpose,
			    &local_j_start_after_transpose,
			    &total_local_size);
    
    if (CPU_Number%2 == 0) {
      if ( (CPU_Rank == CPU_Number/2-1) && (local_i_start+local_Nx-1 < GLOBALNRAD-1) )
	erreur ("ERROR: Bad choice of number of processes for current value of GLOBALNRAD.");
    }
    if ( local_Nx < 2*CPUOVERLAP ) {
      printf ("CPU %d: local_Nx = %d\n", CPU_Rank, local_Nx);
      erreur ("ERROR: One of the processes has a ring number less than the \
total number of ghosts: I must exit.\n");
    }
    one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
    ifront = local_Nx/2 - 1 + (local_Nx%2 == 0 ? 0 : 1);
    /* Each CPU communicates data with its 'CPU_Friend' from the hydro
       mesh to the fftw mesh, and vice versa */
    CPU_Friend = (CPU_Rank + (CPU_Number + one_if_odd)/2)%(CPU_Number + one_if_odd);
    if ( CPU_Number%2 == 1) {
      CPU_NoFriend = (CPU_Number-1)/2;
      CPU_Highest = CPU_NoFriend;
    }
    else {
      CPU_NoFriend = CPU_Number; //does not exist, on purpose. Do not erase!!
      CPU_Highest = CPU_Number-1;
    }
    /* Each CPU communicates boundaries with CPU_Prev and CPU_Next,
       except CPU 0 and CPU_Highest */
    if ( CPU_Rank < (CPU_Number-one_if_odd)/2 ) {
      if ( CPU_Rank == 0 ) {
	CPU_Prev = CPU_Rank;     
	CPU_Next = CPU_Friend;
      } else {
	CPU_Prev = CPU_Friend-1; 
	CPU_Next = CPU_Friend;
      }
    } else {
      if ( CPU_Rank == CPU_Highest ) {
	if ( CPU_Number%2 == 0 ) 
	  CPU_Prev = CPU_Friend;
	else
	  CPU_Prev = 2*CPU_Rank;
	CPU_Next = CPU_Rank; 
      } else {
	CPU_Prev = CPU_Friend;   
	CPU_Next = CPU_Friend+1; 
      }
    }
    /* The hydro mesh dd can now be defined */
    if ( CPU_Rank < (CPU_Number + one_if_odd)/2 ) {
      IMIN = local_i_start;
      IMAX = local_i_start + (local_Nx + (local_Nx%2 == 0 ? 0 : 1))/2 - 1;
      if ( (CPU_Rank == CPU_Highest) && (IMAX >= GLOBALNRAD) )
	IMAX = GLOBALNRAD - 1;
      if ( (CPU_Number%2 == 0) && (CPU_Rank == CPU_Number/2 - 1) )
	transfer_size = 2*(GLOBALNRAD - local_i_start - ifront - 1 + CPUOVERLAP) * 2*(NSEC/2+1);
      else
	transfer_size = 2*(local_Nx - ifront - 1 + CPUOVERLAP) * 2*(NSEC/2 + 1);
      if ( CPU_Rank != CPU_NoFriend ) {
	ffttohydro_transfer = (real *) malloc(sizeof(real) * transfer_size);
	if ( ffttohydro_transfer == NULL )
	  printf("Memory alloc problem with ffttohydro_transfer tab in split.c ...\n");
      }
      sg_split[0] = IMAX;
      sg_split[1] = local_i_start;
      sg_split[2] = local_Nx;
      sg_split[3] = total_local_size;
      sg_split[4] = transfer_size;
      if ( CPU_Rank != CPU_NoFriend )
	MPI_Send(sg_split, 5, MPI_INT, CPU_Friend, 10, MPI_COMM_WORLD);
    } 
    else {
      MPI_Recv(hydro_split, 5, MPI_INT, CPU_Friend, 10, MPI_COMM_WORLD, &fargostat);
      IMAX_friend = hydro_split[0];
      local_i_start_friend = hydro_split[1];
      local_Nx_friend = hydro_split[2];
      total_local_size_friend = hydro_split[3];
      transfer_size_friend = hydro_split[4];
      
      SGP_buffft_Accr_friend = (real *) malloc(sizeof(real)*total_local_size_friend);
      SGP_buffft_Acct_friend = (real *) malloc(sizeof(real)*total_local_size_friend);
      if ( (SGP_buffft_Accr_friend == NULL) || (SGP_buffft_Acct_friend == NULL) )
	printf("Memory alloc problem with SGP_buffft_Acc_friend in split.c ...\n");
      
      ffttohydro_transfer_friend = (real *) malloc(sizeof(real) * transfer_size_friend);
      if ( ffttohydro_transfer_friend == NULL )
	printf("Memory alloc problem with ffttohydro_transfer_friend tab in split.c ...\n");
      
      IMIN = IMAX_friend + 1;
      IMAX = local_i_start_friend + local_Nx_friend - 1;
      if ( (CPU_Rank == CPU_Highest) && (IMAX >= GLOBALNRAD) )
	IMAX = GLOBALNRAD - 1;
    }
    if (CPU_Rank > 0) IMIN -= CPUOVERLAP;
    if (CPU_Rank != CPU_Highest) IMAX += CPUOVERLAP;
    NRAD = IMAX - IMIN + 1;
    if ( NRAD < 2*CPUOVERLAP ) {
      printf ("CPU %d: local_Nx = %d\n", CPU_Rank, NRAD);
      erreur ("ERROR: One of the processes has a ring number less than the \
total number of ghosts: I must exit.\n");
    }
    Zero_or_active = CPUOVERLAP * (CPU_Rank > 0 ? 1 : 0);
    One_or_active = 1+(CPUOVERLAP-1) * (CPU_Rank > 0 ? 1 : 0);
    Max_or_active = NRAD-CPUOVERLAP*(CPU_Rank != CPU_Highest ? 1 : 0);
    MaxMO_or_active = NRAD-1-(CPUOVERLAP-1)*(CPU_Rank != CPU_Highest ? 1 : 0);
    hydro_totalsize = NRAD * NSEC;
    active_hydro_totalsize = (Max_or_active - Zero_or_active + 1)*NSEC;
    
    if ( (CPU_Number%2 == 0) || (CPU_Rank != CPU_NoFriend) ) {
      if ( CPU_Rank >= (CPU_Number+one_if_odd)/2 ) {
	MPI_Ssend (&active_hydro_totalsize, 1, MPI_INT, CPU_Friend, 20, MPI_COMM_WORLD);
	MPI_Ssend (&Zero_or_active, 1, MPI_INT, CPU_Friend, 22, MPI_COMM_WORLD);
      }
      else {
	MPI_Recv (&active_hydro_totalsize_friend, 1, MPI_INT, CPU_Friend, 20, MPI_COMM_WORLD, &fargostat);
	MPI_Recv (&Zero_or_active_friend, 1, MPI_INT, CPU_Friend, 22, MPI_COMM_WORLD, &fargostat);
	dens_friend = (real *) malloc(sizeof(real) * active_hydro_totalsize_friend);
	if ( dens_friend == NULL )  
	  printf("Memory allocation problem with dens_friend in split.c ...\n");
      }
    }
  }
  /* ----------------- */
  /* Let us play debug */
  /* ----------------- */
  if (debug == YES) {
    if ( !SelfGravity || SGZeroMode )
      masterprint ("%d = %d * %d + %d\n", GLOBALNRAD, CPU_Number, size_low, remainder);
    printf ("IMIN process %d : %d\n", CPU_Rank, IMIN);
    printf ("IMAX process %d : %d\n", CPU_Rank, IMAX);
    printf ("NRAD process %d : %d\n", CPU_Rank, NRAD);
    printf ("Zero_or_active process %d : %d\n", CPU_Rank, Zero_or_active);
    printf ("One_or_active process %d : %d\n", CPU_Rank, One_or_active);
    printf ("Max_or_active process %d : %d\n", CPU_Rank, Max_or_active);
    printf ("MaxMO_or_active process %d : %d\n", CPU_Rank, MaxMO_or_active);
    printf ("GLOB process %d : %d\n", CPU_Rank, GLOBALNRAD);
    printf ("CPUOVERLAP process %d : %d\n", CPU_Rank, CPUOVERLAP);
    if ( SelfGravity && !SGZeroMode ) {
      printf ("LocalNx process %d : %d\n", CPU_Rank, local_Nx);
      printf ("LocalIStart process %d : %d\n", CPU_Rank, local_i_start);
      printf ("total_local_size process %d : %d\n", CPU_Rank, total_local_size);
    }
  }
}
