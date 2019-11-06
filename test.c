
#include <stdio.h>
#include <mpi.h>

int main() {
    int rank, size, false=0;
    int right_nbr, left_nbr;
    MPI_Comm ring_comm;
    MPI_Status status[2];
    int value = 1;

    MPI_Init(NULL, NULL);
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Cart_create( MPI_COMM_WORLD, 1, &size, &false, 1, &ring_comm );/*创建一维网格，false意味着没有周期性（两端外的进程标识为MPI_PROC_NULL），1表示可以重排序，得到的拓扑坐标相邻关系为：MPI_PROC_NULL, 0, 1, ... , size-1 , MPI_PROC_NULL*/
    MPI_Cart_shift( ring_comm, 0, 1, &left_nbr, &right_nbr );/*通过在定义的网格上的平移得到左右侧进程的标识*/
    MPI_Comm_rank( ring_comm, &rank );/*得到在新拓扑中的标识或坐标*/
    MPI_Comm_size( ring_comm, &size );/*得到新通信域中进程的个数*/
    MPI_Request req[2];

    if (rank == 0) {/*进程0负责读入数据并向下一个进程传递数据*/
        // scanf( "%d", &value );
        MPI_Isend( &value, 1, MPI_INT, right_nbr, 1, ring_comm, &req[0] );/*将数据传送到右面的进程*/
        MPI_Irecv( &value, 1, MPI_INT, right_nbr, 1, ring_comm, &req[1]);
        MPI_Wait(&req[0], &status[0]);
        printf( "Process %d send %d\n", rank, value );
        MPI_Wait(&req[1], &status[1]);
        value += rank;
        printf( "Process %d got %d\n", rank, value );
    } else {
        MPI_Irecv( &value, 1, MPI_INT, left_nbr, 1, ring_comm, &req[1] );/*后面的进程从左边的进程接收数据*/
        MPI_Isend( &value, 1, MPI_INT, right_nbr, 1, ring_comm, &req[0]);
        MPI_Wait(&req[1], &status[1]);
        value += rank;
        printf( "Process %d got %d\n", rank, value );
        MPI_Wait(&req[0], &status[0]);
        printf( "Process %d send %d\n", rank, value );
        
        

        
        
        
        
        
        
        
        // MPI_Wait(req, status);
        
        
        /*将接收到的数据在传递给右面的进程*/
        
    }
    printf( "Process %d got %d\n", rank, value );/*各进程打印各自得到的数据*/
    // printf("rank%d left is %d, right is %d", rank, left_nbr, right_nbr);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize( );
    return 0;
}
