/*
虚拟拓扑：
	1. 主要的拓扑结构：笛卡尔和图
	2. 虚拟拓扑建立在组和通信因子之上

	MPI_Cart_coords (comm,rank,maxdims,*coords[]) ： 笛卡尔转换函数，把进程秩转换成坐标
	MPI_Cart_create (comm_old,ndims,*dims[],*periods, reorder,*comm_cart) ： 根据笛卡尔拓扑信息创建一个新的通信因子
	MPI_Cart_get (comm,maxdims,*dims,*periods,*coords[]) ： 返回与通信因子相关的笛卡尔拓扑信息
	MPI_Cart_map (comm_old,ndims,*dims[],*periods[],*newrank) ： 返回由MPI系统推荐的笛卡尔映像图
	MPI_Cart_rank (comm,*coords[],*rank) ： 将逻辑进程坐标转换成进程秩
	MPI_Cart_shift (comm,direction,displ,*source,*dest) ： 完成笛卡尔移位操作，根据移位坐标和移位步骤大小确定移位操作
	MPI_Cart_sub (comm,*remain_dims[],*comm_new) ： 该函数用来将通信因子组分成较低维笛卡尔子网格的子组， 并可为每个子组构造一个具 有相关子网格笛卡尔拓扑的通信因子
	MPI_Cartdim_get (comm,*ndims) ： 返回与通信因子相关的笛卡尔结构的维数
	MPI_Dims_create (nnodes,ndims,*dims[]) ： 该函数帮助用户在每个坐标方向选择一种进程的均衡分布， 这取决于被均衡的组中的进 程 数 以 及 由 用 户 指 定 的 限 制 条 件 。 使 用 该 函 数 的 一 种 可 能 方 法 是 将 所 有 进 程 （MPI_COMM_WORLD的组大小）分成一个 n 维的拓扑
	MPI_Graph_create (comm_old,nnodes,*index[],*edges[], reorder,*comm_graph) ：返回一个与图拓扑信息相连的新通信因子
	MPI_Graph_get (comm,maxindex,maxedges,*index[],*edges[]) ： 返回与通信因子相关联的图拓扑信息
	MPI_Graph_map (comm_old,nnodes,*index[],*edges[],*newrank) ： 返回由图的拓扑信息印象
	MPI_Graph_neighbors (comm,rank,maxneighbors,*neighbors[]) ： 返回与指定进程相关的临近节点数组
	MPI_Graphdims_get (comm,*nnodes,*nedges) ： 返回图中的节点数和边数
	MPI_Topo_test (comm,*top_type) ： 返回赋与通信因子的拓扑类型，MPI_GRAPH 图拓扑，MPI_CART 笛卡尔拓扑，MPI_UNDEFINED 没有拓扑



*/
/*
		秩的矩阵
		| 0 | 1 | 2 | 3 |
		| 4 | 5 | 6 | 7 |
		| 8 | 9 | 10| 11|
		| 12| 13| 14| 15|

		秩的坐标及邻近的秩：
		| 秩 | 坐标 | 上 | 下 | 左 | 右 |
		| 0  |(0,0) | -1 | 4  | -1 | 1  |
		| 1  |(0,1) | -1 | 5  | 0  | 2  |
		| 2  |(0,2) | -1 | 6  | 1  | 3  |
		| 3  |(0,3) | -1 | 7  | 2  | -1 |
		| 4  |(1,0) | 0  | 8  | -1 | 5  |
		| 5  |(1,1) | 1  | 9  | 4  | 6  |
		| 6  |(1,2) | 2  | 10 | 5  | 7  |
		| 7  |(1,3) | 3  | 11 | 6  | -1 |
		| 8  |(2,0) | 4  | 12 | -1 | 9  |
		| 9  |(2,1) | 5  | 13 | 8  | 10 |
		| 10 |(2,2) | 6  | 14 | 9  | 11 |
		| 11 |(2,3) | 7  | 15 | 10 | -1 |
		| 12 |(3,0) | 8  | -1 | -1 | 13 |
		| 13 |(3,1) | 9  | -1 | 12 | 14 |
		| 14 |(3,2) | 10 | -1 | 13 | 15 |
		| 15 |(3,3) | 11 | -1 | 14 | -1 |


*/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 
#include "arralloc.h"


#define SIZE 16 
#define UP 0 
#define DOWN 1 
#define LEFT 2 
#define RIGHT 3

int main(int argc, char *argv[])
{
	int numtasks, rank, source, dest, outbuf, i, tag = 1;
	
	int dims[2] = { 4,4 };
	int periods[2] = { 0,0 };
	int reorder = 0;
	int *inbuf;
	inbuf = arralloc(sizeof(int), 1, 4);
	int inf = 0;
	for (inf = 0; inf < 4; inf++) {
		inbuf[inf] = -1;
	}
	// int *nbrs;
	// nbrs = arralloc(sizeof(int), 1, 4);
	int nbrs[4];
	// int *coords;
	// coords = arralloc(sizeof(int), 1, 2);
	int coords[2];
	int id;

	MPI_Request ra, rb, rc, rd, rs;
	// reqs = arralloc(sizeof(MPI_Request), 1, 8);
	MPI_Status stats;
	MPI_Comm cartcomm;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	// int	inbuf[4] = { MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL, };
	
	

	if (numtasks == SIZE) {

		MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);	// 为笛卡尔拓扑创建一个新的通信子

		MPI_Comm_rank(cartcomm, &rank);	// 获取当前进程在通信子中的秩

		MPI_Cart_coords(cartcomm, rank, 2, coords);	// 笛卡尔转换函数，把秩转换成坐标

		//cout << "rank " << rank << " coords : " << coords[0] << "," << coords[1] << endl;

		MPI_Cart_shift(cartcomm, 0, 1, &nbrs[0], &nbrs[1]);	// 笛卡尔移位操作,获取当前坐标位置的上下位置秩

		//cout << "rank " << rank << " vertical shift coords : " << coords[0] << "," << coords[1] << endl;

		MPI_Cart_shift(cartcomm, 1, 1, &nbrs[2], &nbrs[3]);	// 笛卡尔移位操作,获取当前坐标位置的左右位置秩
		MPI_Cart_rank(cartcomm, coords, &id);

		//cout << "rank " << rank << " horizen shift coords : " << coords[0] << "," << coords[1] << endl;
		printf("rank %d coords : %d,%d nbrs: %d, %d, %d, %d\n", id, coords[0], coords[1], nbrs[0], nbrs[1], nbrs[2], nbrs[3]);
		// cout << "rank " << rank << " coords : " << coords[0] << "," << coords[1] << " nbrs : " << nbrs[UP] << "," << nbrs[DOWN] << "," << nbrs[LEFT] << "," << nbrs[RIGHT] << endl;

		outbuf = id;
		// for (i = 0; i < 4; i++) {

			// 将数据分别发送给上下左右4个邻居
			// 从上下左右4个邻居接收数据

			// dest = nbrs[i];
			// source = nbrs[i];

			// if (-1 == dest) continue;

			// 发送
			MPI_Isend(&outbuf, 1, MPI_INT, nbrs[0], tag, MPI_COMM_WORLD, &rs);
			MPI_Isend(&outbuf, 1, MPI_INT, nbrs[1], tag, MPI_COMM_WORLD, &rs);
			MPI_Isend(&outbuf, 1, MPI_INT, nbrs[2], tag, MPI_COMM_WORLD, &rs);
			MPI_Isend(&outbuf, 1, MPI_INT, nbrs[3], tag, MPI_COMM_WORLD, &rs);
			

			// 接收
			MPI_Irecv(&inbuf[0], 1, MPI_INT, nbrs[0], tag, MPI_COMM_WORLD, &ra);
			MPI_Irecv(&inbuf[1], 1, MPI_INT, nbrs[1], tag, MPI_COMM_WORLD, &rb);
			MPI_Irecv(&inbuf[2], 1, MPI_INT, nbrs[2], tag, MPI_COMM_WORLD, &rc);
			MPI_Irecv(&inbuf[3], 1, MPI_INT, nbrs[3], tag, MPI_COMM_WORLD, &rd);
			MPI_Wait(&ra, &stats);
			MPI_Wait(&rb, &stats);
			MPI_Wait(&rc, &stats);
			MPI_Wait(&rd, &stats);
				// printf("rec: %d\n", reqs[i+4]);

		// }
		// printf("%d success!\n", rank_now);
		// MPI_Waitall(8, reqs, &stats);
		// if (rank_now == 0){
			printf("rank= %d coords= %d %d neighbors(u,d,l,r)= %d %d %d %d\n", id, coords[0], coords[1], nbrs[0], nbrs[1], nbrs[2], nbrs[3]);
			printf("rank= %d inbuf(u,d,l,r)= %d %d %d %d\n", id, inbuf[0], inbuf[1], inbuf[2], inbuf[3]);
		// 		}
		

	}
	else
		printf("Must specify %d processors. Terminating.\n", SIZE);
	// int ii;
	// for(ii = 0; ii < SIZE; ii++) {
	// 	if (rank == ii) {
			
	// 	}
	// }

	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
