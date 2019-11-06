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

  int i, j, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc;
  double r;


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

  printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);

  rinit(seed);


  MPI_Comm comm;
  MPI_Status status;

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request reqa, reqb, reqc, reqd;
  MPI_Status sta;
  int dims[1] = {0};
  MPI_Dims_create(size, 1, dims);
  int periods[1] = {1}; 
  int left_nbr, right_nbr;
  MPI_Cart_create( MPI_COMM_WORLD, 1, dims, periods, 0, &comm );
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_shift(comm, 0, 1, &left_nbr, &right_nbr );
  //check size=P
  if (size != 4) {
    printf("Please run with 4 processes.\n");fflush(stdout);
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

  

  MPI_Scatter(&(map[0][0]), M*N, MPI_INT, &(smallmap[0][0]), M*N, MPI_INT, 0, comm);

  


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

  // if (rank == 0) {
  //   for (i = 1; i <= M+2; i++) {
  //     for (j = 1; j <= N+2; j++) {
  //       printf("%3d ", old[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);
  // if (rank == 1) {
  //   for (i = 1; i <= M+2; i++) {
  //     for (j = 1; j <= N+2; j++) {
  //       printf("%3d ", old[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);
  // if (rank == 2) {
  //   for (i = 1; i <= M+2; i++) {
  //     for (j = 1; j <= N+2; j++) {
  //       printf("%3d ", old[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);
  // if (rank == 3) {
  //   for (i = 1; i <= M+2; i++) {
  //     for (j = 1; j <= N+2; j++) {
  //       printf("%3d ", old[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);

  while (step <= maxstep)
    {
      MPI_Isend(&old[M], N, MPI_INT, right_nbr, N, comm, &reqa);
      MPI_Irecv(&old[0], N, MPI_INT, left_nbr, N, comm, &reqb);
      MPI_Isend(&old[1], N, MPI_INT, left_nbr, N, comm, &reqc);
      MPI_Irecv(&old[M+1], N, MPI_INT, right_nbr, N, comm, &reqd);
      MPI_Wait(&reqa, &sta);
      MPI_Wait(&reqb, &sta);
      MPI_Wait(&reqc, &sta);
      MPI_Wait(&reqd, &sta);
      nchange = 0;

      // if (step == 1) {
      //   int ranki;
      //   printf("This is step %d:\n", step);
      //   for (ranki = 0; ranki < 4; ranki++) {
      //     if (rank == ranki) {
      //     for (i = 1; i <= M+2; i++) {
      //       for (j = 1; j <= N+2; j++) {
      //         printf("%3d ", old[i-1][j-1]);
      //       }
      //       printf("\n");
      //     }
      //   }
      //   printf("\n");
      //   MPI_Barrier(comm);
      //   }
      // }

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

  // if (rank == 0) {
  //   for (i = 1; i <= M; i++) {
  //     for (j = 1; j <= N; j++) {
  //       printf("%d ", smallmap[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);
  // if (rank == 1) {
  //   for (i = 1; i <= M; i++) {
  //     for (j = 1; j <= N; j++) {
  //       printf("%d ", smallmap[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);
  // if (rank == 2) {
  //   for (i = 1; i <= M; i++) {
  //     for (j = 1; j <= N; j++) {
  //       printf("%d ", smallmap[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);
  // if (rank == 3) {
  //   for (i = 1; i <= M; i++) {
  //     for (j = 1; j <= N; j++) {
  //       printf("%d ", smallmap[i-1][j-1]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  // MPI_Barrier(comm);

  MPI_Gather(&(smallmap[0][0]), M*N, MPI_INT, &(map[0][0]), M*N, MPI_INT, 0, comm);

  // if (rank == 0){for (i = 0; i < L; i++) {
  //   for (j = 0; j < L; j++) {
  //     printf("%d ", map[i][j]);
  //   }
  //   printf("\n");
  // }}

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

    percwrite("mapdebug.pgm", &(map[0][0]), 8);
  }

  // // MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  free((void *) smallmap);
  return 0;
}
