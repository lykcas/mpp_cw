#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"
#include "arralloc.h"

/*
 * Serial program to test for percolation of a cluster.
 */

int main(int argc, char *argv[])
{
  /*
   *  Define the main arrays for the simulation
   */
  
  
  int old[M+2][N+2], new[M+2][N+2];
  int rank, size;
  int **smallmap;
  smallmap = arralloc(sizeof(int), 2, M, N);
  int **submap;
  submap = arralloc(sizeof(int), 2, L, L);
  int car = L/ (M); //2 * 2
  
  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */

  int **map;
  map = arralloc(sizeof(int), 2, L, L);
  // int smallmap[M][N];
  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, v, o, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc;
  double r;
  
  MPI_Datatype mpi_vec_type;
  MPI_Datatype mpi_rec_vec;
  
  // printf("START2\n");
  int dims[2] = {car, car};
  int periods[2] = {1, 0}; 
  int left_nbr, right_nbr, up_nbr, down_nbr;
  int coord[2];

  
  if (argc != 2)
    {
      printf("Usage: percolate <seed>\n");
      return 1;
    }
  /*
   *  Set most important value: the rock density rho (between 0 and 1)
   */

  rho = 0.411;

  /*
   *  Set the randum number seed and initialise the generator
   */

  seed = atoi(argv[1]);

  if (rank == 0) {printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);}

  rinit(seed);
  MPI_Comm comm;
  // MPI_Status status;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request reqa, reqb, reqc, reqd, reqs;
  MPI_Status sta;
  MPI_Type_vector(M+2, 1, M+2, MPI_INT, &mpi_vec_type);
  MPI_Type_commit(&mpi_vec_type);
  MPI_Type_vector(M, M, L, MPI_INT, &mpi_rec_vec);
  MPI_Type_commit(&mpi_rec_vec);
  MPI_Dims_create(size, 2, dims);
  MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, 0, &comm );
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 2, coord);
  MPI_Cart_shift(comm, 0, 1, &up_nbr, &down_nbr );
  MPI_Cart_shift(comm, 1, 1, &left_nbr, &right_nbr);

  printf("rank %d coords : %d,%d nbrs: up %d, down %d, left %d, right %d\n", rank, coord[0], 
    coord[1], up_nbr, down_nbr, left_nbr, right_nbr);
  if (size != car*car) {
    printf("Please run with 16 processes.\n");fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }


  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */

  if(rank == 0) {
    nhole = 0;

    for (i=0; i < L; i++) {
      for (j=0; j < L; j++) {
        r=uni();
      
        if(r < rho) {
          map[i][j] = 0;
        } else {
          nhole++;
          map[i][j] = nhole;
        }
      }
    }
    printf("percolate: rho = %f, actual density = %f\n",
    rho, 1.0 - ((double) nhole)/((double) L*L) );
  }


  MPI_Bcast(&(map[0][0]), L*L, MPI_INT, 0, comm);


  int x_row = rank / car * (M);
  int x_col = rank % car * (M);

  int i_small = 0, j_small = 0;
  for (i = x_row; i < x_row + (M); i++) {
    for (j = x_col; j < x_col + (M); j++) {
      smallmap[i_small][j_small] = map[i][j];
      j_small += 1;
    }
    i_small += 1;
    j_small = 0;
  }
  


  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */
  for (i=1; i <= M; i++) {
      for (j=1; j <= N; j++) {
	      old[i][j] = smallmap[i-1][j-1];
	    }
    }

   for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][N+1] = 0;
    }

   for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      old[0][j]   = 0;
      old[M+1][j] = 0;
    }


   /*
    *  Update for a fixed number of iterations
    */


  maxstep = 16*L;
  printfreq = 100;

  step = 1;
  nchange = 1;

  while (step <= maxstep)
    {
      
      MPI_Isend(&old[1], (M)+2, MPI_INT, up_nbr, 1, comm, &reqs);
      MPI_Isend(&old[M], (M)+2, MPI_INT, down_nbr, 2, comm, &reqs);
      MPI_Isend(&old[0][1], 1, mpi_vec_type, left_nbr, 1, comm, &reqs);
      MPI_Isend(&old[0][M], 1, mpi_vec_type, right_nbr, 2, comm, &reqs);

      MPI_Irecv(&old[0], (M)+2, MPI_INT, up_nbr, 2, comm, &reqa);
      MPI_Irecv(&old[M+1], (M)+2, MPI_INT, down_nbr, 1, comm, &reqb);
      MPI_Irecv(&old[0][0], 1, mpi_vec_type, left_nbr, 2, comm, &reqc);
      MPI_Irecv(&old[0][M+1], 1, mpi_vec_type, right_nbr, 1, comm, &reqd);
      
      MPI_Wait(&reqa, &sta);
      MPI_Wait(&reqb, &sta);
      MPI_Wait(&reqc, &sta);
      MPI_Wait(&reqd, &sta);
      nchange = 0;

      for (i=1; i<=M; i++) {
	      for (j=1; j<=N; j++) {
          oldval = old[i][j];
          newval = oldval;

          /*
          * Set new[i][j] to be the maximum value of old[i][j]
          * and its four nearest neighbours
          */

          if (oldval != 0) {
            if (old[i-1][j] > newval) newval = old[i-1][j];
            if (old[i+1][j] > newval) newval = old[i+1][j];
            if (old[i][j-1] > newval) newval = old[i][j-1];
            if (old[i][j+1] > newval) newval = old[i][j+1];

            if (newval != oldval)
              {
                ++nchange;
              }
          }

          new[i][j] = newval;
	      }
	    }

      /*
       *  Report progress every now and then
       */

      if (step % printfreq == 0) {
        printf("percolate: number of changes on step %d is %d\n", step, nchange);
	    }

      /*
       *  Copy back in preparation for next step, omitting halos
       */

      for (i=1; i<=M; i++) {
        for (j=1; j<=N; j++)
          {
            old[i][j] = new[i][j];
          }
	    }

      
      step++;

    }

  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. maxstep
   *  is too small)
   */

  if (nchange != 0)
    {
      printf("percolate: WARNING max steps = %d reached before nchange = 0\n", maxstep);
    }

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  
  for (i=1; i<=M; i++) {
    for (j=1; j<=N; j++) {
	    smallmap[i-1][j-1] = old[i][j];
      
	  }
    
  }
  // MPI_Barrier(comm);

  MPI_Isend(&(smallmap[0][0]), (M)*(M), MPI_INT, 0, rank, comm, &reqs);

  
  if (rank == 0) { 
    int rankid = 0;
    for (i = 0; i < L; ) {
      for (j = 0; j < L; ) {
        MPI_Irecv(&(map[i][j]), 1, mpi_rec_vec, MPI_ANY_SOURCE, rankid, comm, &reqa);
        MPI_Wait(&reqa, &sta);
        rankid += 1;
        j += (M);
      }
      i += (M);
    }

    // printf("This is final map:\n");
    // for (i = 0; i < L; i++) {
    //   for (j = 0; j < L; j++) {
    //     printf("%d ", map[i][j]);
    //   }
    //   printf("\n");
    // }
  }




  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */


  if(rank == 0) {
    perc = 0;

    for (itop=0; itop < L; itop++) {
      if (map[itop][L-1] > 0) {
        for (ibot=0; ibot < L; ibot++) {
          if (map[itop][L-1] == map[ibot][0]) {
            perc = 1;
          }
        }
      }
    }

    if (perc != 0)
      {
        printf("percolate: cluster DOES percolate\n");
      }
    else
      {
        printf("percolate: cluster DOES NOT percolate\n");
      }

    /*
    *  Write the map to the file "map.pgm", displaying only the very
    *  largest cluster (or multiple clusters if exactly the same size).
    *  If the last argument here was 2, it would display the largest 2
    *  clusters etc.
    */

    percwrite("map288_16.pgm", &(map[0][0]), 8);
  }

  

  MPI_Finalize();
  free((void *) smallmap);
  free((void *) submap);
  free((void *) map);
  return 0;
}
