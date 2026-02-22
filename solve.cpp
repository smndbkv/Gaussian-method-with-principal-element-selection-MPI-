#include "mpi.h"
#include <time.h>
#include "solve.h"
#include "f.h"

void solve(int argc, char **argv, MPI_Comm com, int p, int k)
{
    int n = 0, m = 0, r = 0, s = 0;
    char *file_name = nullptr;
    double t1 = 0, t2 = 0;

    //    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    if (!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && n > 0 && sscanf(argv[2], "%d", &m) == 1 && m > 0 && sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1 && s >= 0 && s <= 4))
    {
        if (k == 0)
            printf("Usage: %s n m r s <filename>\n", argv[0]);
        return;
    }
    if (argc == 6)
    {
        file_name = argv[5];
    }
    else if (argc == 5 && s == 0)
    {
        if (k == 0)
            printf("Usage: %s n m r 0 filename\n", argv[0]);
        return;
    }
    if (m > n)
    {
        m = n;
    }
    int max_b = (n + m - 1) / m;
    if (p > max_b)
    {
        if (k == 0)
            printf("The number of processes exceeds the number of block rows\n");
        return;
    }
    // int rows = get_rows(n, m, p, k);
    int max_b_rows = get_max_block_rows(n, m, p, k);

    // printf("rows in process %d = %d\n", k, rows);
    double *a = new double[max_b_rows * m * n],
           *buf = new double[m * n];
    if (argc == 6)
    {
        int err = read_matrix(a, n, m, p, k, file_name, buf, com);
        if (err)
        {
            if (k == 0)
                printf("Cannot read file %s\n", file_name);
            delete[] a;
            delete[] buf;
            return;
        }
    }
    else
    {
        init_matrix(a, n, m, p, k, s);
    }
    if (k == 0)
        printf("Initial matrix:\n");
    print_matrix(a, n, m, p, k, buf, r, com);

    // printf("max_block_rows = %d\n", max_b_rows);
    double *b = new double[max_b_rows * m], *x_exact = new double[max_b_rows * m];
    init_x_exact(x_exact, n, m, p, k);
    init_b(a, b, n, m, p, k);

    if (k == 0)
        printf("Initial exact solution:\n");
    print_vector(x_exact, n, m, p, k, buf, r, com);
    if (k == 0)
        printf("Initial right side:\n");
    print_vector(b, n, m, p, k, buf, r, com);

    double *x = new double[max_b_rows * m];

    t2 = clock();
    gaussian_method(a, b, x, n, m, p, k, buf, com);
    t2 = (clock() - t2) / CLOCKS_PER_SEC;

    t2 = clock();
    double r1 = get_r1(a, x_exact /* = x*/, b, buf, n, m, p, k, com), r2 = get_r2(x_exact /* = x*/, x_exact, n, m, p, k, com);
    t2 = (clock() - t2) / CLOCKS_PER_SEC;

    if (k == 0)
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
               argv[0], 11, r1, r2, t1, t2, s, n, m, p);

    delete[] a;
    delete[] buf;
    delete[] b;
    delete[] x_exact;
}