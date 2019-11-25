#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>



void broadcast(int **map, int **smallmap, int M, int N, int l, int rank,
               int dims[2]) {
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
               int down_nbr, int left_nbr, int right_nbr) {
  int maxstep = 16 * l;
  int printfreq = 100;
  int step = 1;
  int nchange = 1, nchange_count = 0;
  int oldval, newval, i, j, nchange_sum;

  MPI_Request reqs, reqa, reqb, reqc, reqd;
  MPI_Status sta;
  MPI_Datatype mpi_vec_type;
  MPI_Type_vector(M + 2, 1, N + 2, MPI_INT, &mpi_vec_type);
  MPI_Type_commit(&mpi_vec_type);

  while (step <= maxstep) {
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

    /* Calculate the sum of the numbers in map. */
    MPI_Reduce(&local_mapsum, &global_mapsum, 1, MPI_LONG_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);

    /* Calculate the average of the number in map. */
    
    /* Find whether the processors are changing or stopping changing. */
    nchange_count = nchange;
    MPI_Allreduce(&nchange_count, &nchange_sum, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    /* If the processors all stop changing, stop the iteration part. */
    if (nchange_sum == 0) {
      if (rank == 0) printf("%d\n", step);
      break;
    }

    if (step % printfreq == 0) {
      printf("percolate: number of changes on step %d is %d\n", step, nchange);
      if (rank == 0) {
        printf("The average of the map array on step %d id %.4lf\n", step,
              (double)global_mapsum / (double)(l * l));
      }
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