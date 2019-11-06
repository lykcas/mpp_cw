#include <stdio.h>
#include <mpi.h>

int main() {
    int rank, size, i;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int data[4] = {1, 2, 3, 4};
    int value;
    MPI_Status status;
    //MPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //data += rank;
    if (rank == 0) {
        MPI_Scatter(data, 1, MPI_INT, &value, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } 


    printf("data = %d in %d process\n", value, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}