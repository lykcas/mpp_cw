

void broadcast(int **map, int **smallmap, int M, int N, int l, int rank, int dims[2]);

void halo(int **old, int **smallmap, int M, int N);

void iteration(int **old, int **new, int M, int N, int l, int rank, int up_nbr, int down_nbr,
               int left_nbr, int right_nbr);

void gather(int **smallmap, int **old, int **map, int M, int N, int l,
            int rank);

