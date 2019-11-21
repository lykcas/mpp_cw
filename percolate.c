#include "percolate.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "arralloc.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: percolate <seed>\n");
    return 1;
  }

  int seed;
  double rho;
  int i, j, v, o, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc, nchange_count, nchange_sum;
  double r;
  long long local_mapsum, global_mapsum;

  int dims[2] = {0, 0};
  int periods[2] = {1, 0};
  int left_nbr, right_nbr, up_nbr, down_nbr;
  int coord[2];
  struct time *time_info;
  time_info = (struct time *)arralloc(sizeof(struct time), 1, 1);
  // double start_global, start_iteration, end_global, end_iteration;
  // double start_commu, end_commu, time_commu = 0, start_update, end_update,
  //                                time_update = 0;
  int rank, size;

  rho = 0.411;

  seed = atoi(argv[1]);

  if (rank == 0) {
    printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);
  }

  rinit(seed);

  MPI_Comm comm;

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // MPI_Request reqa, reqb, reqc, reqd, reqs;
  // MPI_Status sta;
  MPI_Dims_create(size, 2, dims);
  int M = L / dims[0], N = L / dims[1];
  MPI_Datatype mpi_vec_type;
  MPI_Type_vector(M + 2, 1, N + 2, MPI_INT, &mpi_vec_type);
  MPI_Type_commit(&mpi_vec_type);
  MPI_Datatype mpi_rec_vec;
  MPI_Type_vector(M, N, L, MPI_INT, &mpi_rec_vec);
  MPI_Type_commit(&mpi_rec_vec);
  // int old[M + 2][N + 2], new[M + 2][N + 2];
  int **old;
  old = arralloc(sizeof(int), 2, M + 2, N + 2);
  int **new;
  new = arralloc(sizeof(int), 2, M + 2, N + 2);
  int **smallmap;
  smallmap = arralloc(sizeof(int), 2, M, N);
  int **submap;
  submap = arralloc(sizeof(int), 2, L, L);
  int **map;
  map = arralloc(sizeof(int), 2, L, L);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 2, coord);
  MPI_Cart_shift(comm, 0, 1, &up_nbr, &down_nbr);
  MPI_Cart_shift(comm, 1, 1, &left_nbr, &right_nbr);

  // printf("rank %d coords : %d,%d nbrs: up %d, down %d, left %d, right %d\n",
  //        rank, coord[0], coord[1], up_nbr, down_nbr, left_nbr, right_nbr);
  if (size != dims[0] * dims[1]) {
    printf("Please run with %d processes.\n", dims[0] * dims[1]);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  time_info->start_global = MPI_Wtime();
  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */
  initialmap((int **)map, rank, L, rho);

  broadcast((int **)map, (int **)smallmap, M, N, L, rank, dims);

  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */
  halo((int **)old, (int **)smallmap, M, N);

  /*
   *  Update for a fixed number of iterations
   */
  iteration((int **)old, (int **)new, M, N, L, rank, up_nbr,
            down_nbr, left_nbr, right_nbr, time_info);

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  gather((int **)smallmap, (int **)old, (int **)map, M, N, L, rank);

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */
  judge_write((int **)map, rank, L);

  time_info->end_global = MPI_Wtime();
  MPI_Finalize();

  time_print(rank, time_info);

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

void broadcast(int **map, int **smallmap, int M, int N, int l, int rank, int dims[2]) {
  MPI_Bcast(&(map[0][0]), l * l, MPI_INT, 0, MPI_COMM_WORLD);
  int x_row = (rank) / dims[1] * (M);
  int x_col = (rank) % dims[1] * (N);
  int i_small = 0, j_small = 0, i, j;
  for (i = x_row; i < x_row + (M); i++) {
    for (j = x_col; j < x_col + (N); j++) {
      smallmap[i_small][j_small] = map[i][j];
      j_small += 1;
    }
    i_small += 1;
    j_small = 0;
  }
}

void halo(int **old, int **smallmap, int M, int N) {
  int i, j;
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= N; j++) {
      old[i][j] = smallmap[i - 1][j - 1];
    }
  }

  for (i = 0; i <= M + 1; i++) {
    old[i][0] = 0;
    old[i][N + 1] = 0;
  }

  for (j = 0; j <= N + 1; j++) {
    old[0][j] = 0;
    old[M + 1][j] = 0;
  }
}

void iteration(int **old, int **new, int M, int N, int l, int rank, int up_nbr,
               int down_nbr, int left_nbr, int right_nbr,
               struct time *time_info) {
  int maxstep = 16 * L;
  int printfreq = 50;
  int step = 1;
  int nchange = 1, nchange_count = 0;
  int oldval, newval, i, j, nchange_sum;
  double temp_commu = 0, temp_update = 0;
  
  MPI_Request reqs, reqa, reqb, reqc, reqd;
  MPI_Status sta;
  MPI_Datatype mpi_vec_type;
  MPI_Type_vector(M + 2, 1, N + 2, MPI_INT, &mpi_vec_type);
  MPI_Type_commit(&mpi_vec_type);
  time_info->start_iteration = MPI_Wtime();
  while (step <= maxstep) {
    double start_commu = MPI_Wtime();
    MPI_Isend(&(old[1][0]), (N) + 2, MPI_INT, up_nbr, 1, MPI_COMM_WORLD, &reqs);
    MPI_Isend(&(old[M][0]), (N) + 2, MPI_INT, down_nbr, 2, MPI_COMM_WORLD,
              &reqs);
    MPI_Isend(&(old[0][1]), 1, mpi_vec_type, left_nbr, 1, MPI_COMM_WORLD,
              &reqs);
    MPI_Isend(&(old[0][N]), 1, mpi_vec_type, right_nbr, 2, MPI_COMM_WORLD,
              &reqs);

    MPI_Irecv(&(old[0][0]), (N) + 2, MPI_INT, up_nbr, 2, MPI_COMM_WORLD, &reqa);
    MPI_Irecv(&(old[M + 1][0]), (N) + 2, MPI_INT, down_nbr, 1, MPI_COMM_WORLD,
              &reqb);
    MPI_Irecv(&(old[0][0]), 1, mpi_vec_type, left_nbr, 2, MPI_COMM_WORLD,
              &reqc);
    MPI_Irecv(&(old[0][N + 1]), 1, mpi_vec_type, right_nbr, 1, MPI_COMM_WORLD,
              &reqd);

    MPI_Wait(&reqa, &sta);
    MPI_Wait(&reqb, &sta);
    MPI_Wait(&reqc, &sta);
    MPI_Wait(&reqd, &sta);
    double end_commu = MPI_Wtime();
    temp_commu += (end_commu - start_commu);
    double start_update = MPI_Wtime();
    nchange = 0;
    long long local_mapsum = 0;
    long long global_mapsum = 0;
    for (i = 1; i <= M; i++) {
      for (j = 1; j <= N; j++) {
        oldval = old[i][j];
        newval = oldval;
        local_mapsum += oldval;

        /*
         * Set new[i][j] to be the maximum value of old[i][j]
         * and its four nearest neighbours
         */

        if (oldval != 0) {
          if (old[i - 1][j] > newval) newval = old[i - 1][j];
          if (old[i + 1][j] > newval) newval = old[i + 1][j];
          if (old[i][j - 1] > newval) newval = old[i][j - 1];
          if (old[i][j + 1] > newval) newval = old[i][j + 1];

          if (newval != oldval) {
            ++nchange;
          }
        }

        new[i][j] = newval;
      }
    }
    double end_update = MPI_Wtime();
    temp_update += (end_update - start_update);
    

    if (step % printfreq == 0) {
      MPI_Reduce(&local_mapsum, &global_mapsum, 1, MPI_LONG_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
      // printf("percolate: number of changes on step %d is %d\n", step, nchange);
      // if (rank == 0)
      //   printf("The average of the map array on step %d id %.4lf\n", step,
      //          (double)global_mapsum / (double)(L * L));
    }
    nchange_count = nchange;
    MPI_Allreduce(&nchange_count, &nchange_sum, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    if (nchange_sum == 0) {
      if (rank == 0)
        printf("Percolate finish. There are no changes in the map array.\n");
      break;
    }

    /*
     *  Copy back in preparation for next step, omitting halos
     */

    for (i = 1; i <= M; i++) {
      for (j = 1; j <= N; j++) {
        old[i][j] = new[i][j];
      }
    }

    step++;
  }
  time_info->end_iteration = MPI_Wtime();
  time_info->time_commu = temp_commu;
  time_info->time_update = temp_update;

  if (nchange != 0) {
    printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
           maxstep);
  }
}

void gather(int **smallmap, int **old, int **map, int M, int N, int l,
            int rank) {
  MPI_Request reqs, reqa;
  MPI_Status sta;
  MPI_Datatype mpi_rec_vec;
  MPI_Type_vector(M, N, l, MPI_INT, &mpi_rec_vec);
  MPI_Type_commit(&mpi_rec_vec);
  int i, j;
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= N; j++) {
      smallmap[i - 1][j - 1] = old[i][j];
    }
  }

  MPI_Isend(&(smallmap[0][0]), (M) * (N), MPI_INT, 0, rank, MPI_COMM_WORLD,
            &reqs);

  if (rank == 0) {
    int rankid = 0;
    for (i = 0; i < l;) {
      for (j = 0; j < l;) {
        MPI_Irecv(&(map[i][j]), 1, mpi_rec_vec, MPI_ANY_SOURCE, rankid,
                  MPI_COMM_WORLD, &reqa);
        MPI_Wait(&reqa, &sta);
        rankid += 1;
        j += (N);
      }
      i += (M);
    }
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
    percwrite("map288.pgm", map, 8);
  }
}

void time_print(int rank, struct time *time_info) {
  printf(
      "%d:\nCommunicaiton time is %.4f.\nUpdating time is %.4f.\nIteration "
      "time is %.4f.\nRuntime is %.4f.\n",
      rank, time_info->time_commu, time_info->time_update,
      time_info->end_iteration - time_info->start_iteration,
      time_info->end_global - time_info->start_global);
}
