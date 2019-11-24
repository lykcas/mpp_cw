/*
 *  Main header file for percolation code.
 */

/*
 *  System size L
 */

#define L 1152

/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over P
 *  processes, then M = L/P and N = L
 */
struct time {
  double start_global;
  double end_global;
  double start_iteration;
  double end_iteration;
  double time_commu;
  double time_update;
};

void initialmap(int **map, int rank, int l, double rho);

void broadcast(int **map, int **smallmap, int M, int N, int l, int rank, int dims[2]);

void halo(int **old, int **smallmap, int M, int N);

void iteration(int **old, int **new, int M, int N, int l, int rank, int up_nbr, int down_nbr,
               int left_nbr, int right_nbr, struct time *time_info);

void gather(int **smallmap, int **old, int **map, int M, int N, int l,
            int rank);

void judge_write(int **map, int rank, int l);

void time_print(int rank, struct time *time_info);

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void percwrite(char *percfile, int **map, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
