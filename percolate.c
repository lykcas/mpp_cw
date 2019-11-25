#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "arralloc.h"
#include "percwrite.h"
#include "percolate.h"
#include "grid.h"

int main(int argc, char *argv[]) {

  int seed, len;
  if (argc < 2) {
    printf("Usage: percolate <seed>\n");
    return 1;
  } else if (argc == 2) {
    len = 288;
  } else {
    len = atoi(argv[2]);
  }
  
  double rho;
  int i, j, v, o, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc, nchange_count, nchange_sum;
  double r;
  long long local_mapsum, global_mapsum;

  int dims[2] = {0, 0};
  int periods[2] = {1, 0};
  int left_nbr, right_nbr, up_nbr, down_nbr;
  int coord[2];
  int rank, size;

  rho = 0.411;
  
  seed = atoi(argv[1]);

  if (rank == 0) {
    printf("percolate: params are len = %d, rho = %f, seed = %d\n", len, rho, seed);
  }

  rinit(seed);

  MPI_Comm comm;

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  int M = len / dims[0], N = len / dims[1];
  MPI_Datatype mpi_vec_type;
  MPI_Type_vector(M + 2, 1, N + 2, MPI_INT, &mpi_vec_type);
  MPI_Type_commit(&mpi_vec_type);
  MPI_Datatype mpi_rec_vec;
  MPI_Type_vector(M, N, len, MPI_INT, &mpi_rec_vec);
  MPI_Type_commit(&mpi_rec_vec);

  int **old;
  old = arralloc(sizeof(int), 2, M + 2, N + 2);
  int **new;
  new = arralloc(sizeof(int), 2, M + 2, N + 2);
  int **smallmap;
  smallmap = arralloc(sizeof(int), 2, M, N);
  int **submap;
  submap = arralloc(sizeof(int), 2, len, len);
  int **map;
  map = arralloc(sizeof(int), 2, len, len);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 2, coord);
  MPI_Cart_shift(comm, 0, 1, &up_nbr, &down_nbr);
  MPI_Cart_shift(comm, 1, 1, &left_nbr, &right_nbr);

  if (size != dims[0] * dims[1]) {
    printf("Please run with %d processes.\n", dims[0] * dims[1]);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */
  initialmap((int **)map, rank, len, rho);

  broadcast((int **)map, (int **)smallmap, M, N, len, rank, dims);

  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */
  halo((int **)old, (int **)smallmap, M, N);

  /*
   *  Update for a fixed number of iterations
   */
  double start_test = MPI_Wtime();
  iteration((int **)old, (int **)new, M, N, len, rank, up_nbr, down_nbr, left_nbr,
            right_nbr);
  double end_test = MPI_Wtime();
  if (rank == 0) printf("%f\n", end_test - start_test);

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  gather((int **)smallmap, (int **)old, (int **)map, M, N, len, rank);

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */
  judge_write((int **)map, rank, len);

  MPI_Finalize();


  free((void *)smallmap);
  free((void *)submap);
  free((void *)map);
  free((void *)old);
  free((void *)new);
  return 0;
}

void initialmap(int **map, int rank, int l, double rho) {
  if (rank == 0) {
    int nhole = 0, i, j;
    double r;
    for (i = 0; i < l; i++) {
      for (j = 0; j < l; j++) {
        r = uni();

        if (r < rho) {
          map[i][j] = 0;
        } else {
          nhole++;
          map[i][j] = nhole;
        }
      }
    }
    printf("percolate: rho = %f, actual density = %f\n", rho,
           1.0 - ((double)nhole) / ((double)l * l));
  }
}

void judge_write(int **map, int rank, int l) {
  int perc, itop, ibot;
  if (rank == 0) {
    perc = 0;
    for (itop = 0; itop < l; itop++) {
      if (map[itop][l - 1] > 0) {
        for (ibot = 0; ibot < l; ibot++) {
          if (map[itop][l - 1] == map[ibot][0]) {
            perc = 1;
          }
        }
      }
    }
    if (perc != 0) {
      printf("percolate: cluster DOES percolate\n");
    } else {
      printf("percolate: cluster DOES NOT percolate\n");
    }
    percwrite("map.pgm", map, 8, l);
  }
}

