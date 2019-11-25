/*
 *  Main header file for percolation code.
 */

/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over P
 *  processes, then M = L/P and N = L
 */

void initialmap(int **map, int rank, int l, double rho);

void judge_write(int **map, int rank, int l);

void rinit(int ijkl);

float uni(void);
