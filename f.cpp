#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <sys/sysinfo.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include "io_status.h"
#include "f.h"

double get_full_time()
{
    struct timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + t.tv_usec / 1e6;
}
double get_cpu_time()
{
    struct rusage t;
    getrusage(RUSAGE_THREAD, &t);
    return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1e6;
}

int g2l_b(int, int, int p, int, int i_glob) // block to block
{
    return i_glob / p;
}
int l2g_b(int, int, int p, int k, int i_loc) // block to block
{
    return i_loc * p + k;
}
int l2g(int, int m, int p, int k, int i_loc) // double to double
{
    int i_loc_m = i_loc / m;
    int i_glob_m = i_loc_m * p + k;
    return i_glob_m * m + i_loc % m;
}

int get_block_rows(int n, int m, int p, int k)
{
    int b = (n + m - 1) / m;
    return (k < b % p ? b / p + 1 : b / p);
}

int get_rows(int n, int m, int p, int k)
{
    int b = (n + m - 1) / m;
    int b_last = (b - 1) % p;
    int b_loc = (b % p <= k ? b / p : b / p + 1);
    if (k != b_last)
    {
        return b_loc * m;
    }
    int l = n % m;
    if (l == 0)
    {
        return b_loc * m;
    }
    if (b_loc == 0)
    {
        return 0;
    }

    return (b_loc - 1) * m + l;
}
// int get_rows(int n, int m, int p, int k)
// {
//     int nb = (n + m - 1) / m;
//     int last_block = nb - 1;
//     int last_owner = last_block % p;

//     int count = (nb - k + p - 1) / p;

//     if (count == 0)
//         return 0;

//     int total = count * m;
//     if (last_owner == k)
//     {
//         total = total - m + (n - last_block * m);
//     }

//     return total;
// }

int get_max_block_rows(int n, int m, int p, int)
{
    int b = (n + m - 1) / m;
    return (b + p - 1) / p;
}

double f(int s, int n, int i, int j)
{
    switch (s)
    {
    case 1:
        return n - std::max(i, j) + 1;
        break;
    case 2:
        return std::max(i, j);
        break;
    case 3:
        return std::abs(i - j);
        break;
    case 4:
        return 1. / (i + j - 1);
        break;
    }
    return 0;
}
void init_matrix(double *a, int n, int m, int p, int k, int s)
{
    int i_loc, j_loc, i_glob, j_glob;
    int rows = get_rows(n, m, p, k);
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        i_glob = l2g(n, m, p, k, i_loc);
        for (j_loc = 0; j_loc < n; j_loc++)
        {
            j_glob = j_loc;
            a[i_loc * n + j_loc] = f(s, n, i_glob + 1, j_glob + 1);
        }
    }
}

int print_array(double *a, int n, int m, int printed_rows, int max_print)
{
    if (printed_rows >= max_print)
    {
        return 0;
    }
    int p_n = (n > max_print ? max_print : n);
    int p_m = (printed_rows + m <= max_print ? m : max_print - printed_rows);
    for (int i = 0; i < p_m; i++)
    {
        for (int j = 0; j < p_n; j++)
        {
            printf(" %10.3e", a[i * n + j]);
        }
        printf("\n");
    }
    return p_m;
}

void print_matrix(double *a, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm com)
{
    int main_k = 0; // только процесс с номером ноль имеет достув к stdout
    int max_b = (n + m - 1) / m, printed_rows = 0;
    for (int b = 0; b < max_b; b++)
    {
        int owner = b % p, b_loc = b / p, rows = (m <= n - b * m ? m : n - b * m);
        if (k == main_k)
        {
            if (owner == main_k)
            {
                printed_rows += print_array(a + b_loc * n * m, n, rows, printed_rows, max_print);
            }
            else
            {
                MPI_Status st;
                MPI_Recv(buf, n * rows, MPI_DOUBLE, owner, 0, com, &st);
                printed_rows += print_array(buf, n, rows, printed_rows, max_print);
            }
        }
        else if (owner == k)
        {
            MPI_Send(a + b_loc * n * m, n * rows, MPI_DOUBLE, main_k, 0, com);
        }
    }
}

int print_vector_part(double *v, int rows, int printed, int max_print)
{
    if (printed >= max_print)
        return 0;
    int to_print = rows;
    if (printed + rows > max_print)
        to_print = max_print - printed;
    for (int i = 0; i < to_print; i++)
        printf(" %10.3e", v[i]);
    return to_print;
}

void print_vector(double *b, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm com)
{
    int main_k = 0;
    int max_b = (n + m - 1) / m;
    int printed = 0;
    for (int b_idx = 0; b_idx < max_b; b_idx++)
    {
        int owner = b_idx % p;
        int b_loc = b_idx / p;
        int rows = (b_idx * m + m <= n) ? m : (n - b_idx * m);
        if (k == main_k)
        {
            if (owner == main_k)
            {
                printed += print_vector_part(b + b_loc * m, rows, printed, max_print);
            }
            else
            {
                MPI_Status st;
                MPI_Recv(buf, rows, MPI_DOUBLE, owner, 0, com, &st);
                printed += print_vector_part(buf, rows, printed, max_print);
            }
        }
        else if (owner == k)
        {
            MPI_Send(b + b_loc * m, rows, MPI_DOUBLE, main_k, 0, com);
        }
    }
    if (k == main_k)
        printf("\n");
}

int read_array(FILE *fp, double *a, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (fscanf(fp, "%lf", a + i) != 1)
        {
            return 1;
        }
        // printf("%lf ", a[i]);
    }
    // printf("\n");
    return 0;
}

int read_matrix(double *a, int n, int m, int p, int k, const char *name, double *buf, MPI_Comm com)
{
    int main_k = 0;
    FILE *fp = nullptr;
    int err = 0;
    if (k == main_k)
    {
        fp = fopen(name, "r");
        if (fp == nullptr)
            err = 1;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, com);
    if (err)
        return err;
    memset(buf, 0, n * m * sizeof(double));
    int b, max_b = (n + m - 1) / m;
    for (b = 0; b < max_b; b++)
    {
        int owner = b % p;
        int rows = (b * m + m <= n ? m : n - b * m);
        // int rows = get_rows(n, m, p, k);
        int b_loc = b / p;
        if (k == main_k)
        {
            err += read_array(fp, buf, n * rows);
            if (owner == main_k)
            {
                memcpy(a + b_loc * n * m, buf, n * rows * sizeof(double));
            }
            else
            {
                MPI_Send(buf, n * rows, MPI_DOUBLE, owner, 0, com);
            }
        }
        else if (owner == k)
        {
            MPI_Status st;
            MPI_Recv(a + b_loc * m * n, n * rows, MPI_DOUBLE, main_k, 0, com, &st);
        }
    }

    if (k == main_k)
    {
        fclose(fp);
        fp = nullptr;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, com);
    if (err)
        return err;
    return 0;
}

void matrix_mult_vector(double *a, double *b, double *c, int n, int m, int p, int k, MPI_Comm com)
{
    int b_rows = get_block_rows(n, m, p, k);
    int rows = get_rows(n, m, p, k);
    int max_b_rows = get_max_block_rows(n, m, p, k);
    int dst = (k - 1 + p) % p;
    int src = (k + 1 + p) % p;

    memset(c, 0, rows * sizeof(double));

    for (int s = 0; s < p; s++)
    {
        int sk = (k + s) % p;
        int sk_b_rows = get_block_rows(n, m, p, sk);

        for (int i = 0; i < b_rows; i++)
        {
            int i_glob = l2g_b(n, m, p, k, i);
            int row_h = (i_glob * m + m <= n) ? m : (n - i_glob * m);

            for (int sk_i = 0; sk_i < sk_b_rows; sk_i++)
            {
                int sk_i_glob = l2g_b(n, m, p, sk, sk_i);
                int w = (sk_i_glob * m + m <= n) ? m : (n - sk_i_glob * m);

                for (int ii = 0; ii < row_h; ii++)
                {
                    double sum = 0.0;
                    int local_row = i * m + ii;
                    for (int jj = 0; jj < w; jj++)
                    {
                        int col_glob = sk_i_glob * m + jj;
                        sum += a[local_row * n + col_glob] * b[sk_i * m + jj];
                    }
                    c[local_row] += sum;
                }
            }
        }

        MPI_Status st;
        MPI_Sendrecv_replace(b, max_b_rows * m, MPI_DOUBLE, dst, 0, src, 0, com, &st);
    }
}

void init_b(double *a, double *b, int n, int m, int p, int k)
{
    int i_loc, j_loc;
    int rows = get_rows(n, m, p, k);
    double sum = 0;
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        sum = 0;
        for (j_loc = 0; j_loc <= (n - 1) / 2; j_loc++)
        {
            sum += a[i_loc * n + 2 * j_loc];
        }
        // printf(" k = %d, sum = %lf\n", k, sum);
        b[i_loc] = sum;
    }
}

void init_x_exact(double *x, int n, int m, int p, int k)
{
    int i_loc, i_glob;
    int rows = get_rows(n, m, p, k);
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        i_glob = l2g(n, m, p, k, i_loc);
        x[i_loc] = (i_glob + 1) % 2;
    }
}

double l1_norm(double *x, int n, int m, int p, int k, MPI_Comm com)
{
    int i_loc;
    int rows = get_rows(n, m, p, k);
    // printf("%d\n", rows);
    double sum_loc = 0, sum_glob = 0;
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        sum_loc += fabs(x[i_loc]);
    }
    MPI_Allreduce(&sum_loc, &sum_glob, 1, MPI_DOUBLE, MPI_SUM, com);
    return sum_glob;
}

double l2_norm_matrix(double *a, int n, int m, int p, int k, MPI_Comm com)
{
    double sum = 0, sum_glob;
    int rows = get_rows(n, m, p, k);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum += a[i * n + j] * a[i * n + j];
        }
    }
    MPI_Allreduce(&sum, &sum_glob, 1, MPI_DOUBLE, MPI_SUM, com);

    return sqrt(sum_glob);
}
double get_r1(double *a, double *x, double *b, double *buf, int n, int m, int p, int k, MPI_Comm com)
{
    matrix_mult_vector(a, x, buf, n, m, p, k, com);
    int i_loc;
    int rows = get_rows(n, m, p, k);
    // printf("%d\n", rows);
    double sum_loc1 = 0, sum_glob1 = 0, sum_loc2 = 0, sum_glob2 = 0;
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        sum_loc1 += fabs(buf[i_loc] - b[i_loc]);
        sum_loc2 += fabs(b[i_loc]);
    }
    MPI_Allreduce(&sum_loc1, &sum_glob1, 1, MPI_DOUBLE, MPI_SUM, com);
    MPI_Allreduce(&sum_loc2, &sum_glob2, 1, MPI_DOUBLE, MPI_SUM, com);
    if (CMP(sum_glob2, 0))
        return -1;
    return sum_glob1 / sum_glob2;
}

double get_r2(double *x, double *x_exact, int n, int m, int p, int k, MPI_Comm com)
{
    int i_loc;
    int rows = get_rows(n, m, p, k);
    // printf("%d\n", rows);
    double sum_loc1 = 0, sum_glob1 = 0;
    double sum_loc2 = 0, sum_glob2 = 0;
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        sum_loc1 += fabs(x[i_loc] - x_exact[i_loc]);
        sum_loc2 += fabs(x[i_loc]);
    }
    MPI_Allreduce(&sum_loc1, &sum_glob1, 1, MPI_DOUBLE, MPI_SUM, com);
    MPI_Allreduce(&sum_loc2, &sum_glob2, 1, MPI_DOUBLE, MPI_SUM, com);
    if (CMP(sum_glob2, 0))
    {
        return -1;
    }
    return sum_glob1 / sum_glob2;
}

void get_block(double *a, int n, int m, int p, int k, int i, int j, double *c, int &v, int &h)
{
    int r = 0, t = 0;
    int kk = n / m;
    int l = n - kk * m;
    int i_glob = l2g_b(n, m, p, k, i), j_glob = j;
    v = (i_glob < kk ? m : l);
    h = (j_glob < kk ? m : l);
    // double *pa = a + i * n * m + j * m;
    for (r = 0; r < v; r++)
    {
        for (t = 0; t < h; t++)
        {
            c[r * h + t] = a[(i * m + r) * n + j * m + t];
        }
    }
}

void get_block_buf(double *a, int n, int m, int p, int k, int s_loc, int j, double *c, int &v, int &h)
{
    int r = 0, t = 0;
    int kk = n / m;
    int l = n - kk * m;
    int i_glob = l2g_b(n, m, p, k, s_loc), j_glob = j;
    v = (i_glob < kk ? m : l);
    h = (j_glob < kk ? m : l);
    // double *pa = a + i * n * m + j * m;
    int i = 0;
    for (r = 0; r < v; r++)
    {
        for (t = 0; t < h; t++)
        {
            c[r * h + t] = a[(i * m + r) * n + j * m + t];
        }
    }
}
void set_block(double *a, int n, int m, int i, int j, double *c, int v, int h)
{
    for (int r = 0; r < v; r++)
    {
        for (int t = 0; t < h; t++)
        {
            a[(i * m + r) * n + (j * m + t)] = c[r * h + t];
        }
    }
}

void get_block_vector(double *b, int n, int m, int p, int k, int i, double *c, int &h)
{
    int i_glob = l2g_b(n, m, p, k, i);
    h = (i_glob < n / m ? m : n % m);
    for (int r = 0; r < h; r++)
    {
        c[r] = b[m * i + r];
    }
}
void set_block_vector(double *b, int m, int i, double *c, int h)
{
    for (int r = 0; r < h; r++)
    {
        b[m * i + r] = c[r];
    }
}

void print_matrix_local(double *c, int v, int h, int k, int main_k = 0)
{
    for (int i = 0; i < v; i++)
    {
        for (int j = 0; j < h; j++)
        {
            if (k == main_k)
            {
                printf(" %10.3e", c[i * h + j]);
            }
        }
        if (k == main_k)
            printf("\n");
    }
}

double matrix_norm(double *a, int n, int m, int p, int k, MPI_Comm com)
{
    int i, j;
    double max = -1, sum = 0;
    int rows = get_rows(n, m, p, k);
    for (i = 0; i < rows; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++)
        {
            sum += fabs(a[i * n + j]);
        }
        if (max < sum)
        {
            max = sum;
        }
    }
    double max_glob;
    MPI_Allreduce(&max, &max_glob, 1, MPI_DOUBLE, MPI_MAX, com);
    return max_glob;
}

double matrix_norm_local(double *c, int v)
{
    int i, j;
    double max = -1, sum = 0;
    for (i = 0; i < v; i++)
    {
        sum = 0;
        for (j = 0; j < v; j++)
        {
            sum += fabs(c[i * v + j]);
        }
        if (max < sum)
        {
            max = sum;
        }
    }
    return max;
}

bool inverse(double *a, int n, double *c, double nrm_a)
{
    // const double nrm = norm(a, n);
    const double nrm = nrm_a;
    const double eps = EPS * nrm;
    double max, d;
    int i, j, k;
    int max_i;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            c[i * n + j] = (i == j) ? 1 : 0;
        }
    }

    for (i = 0; i < n; ++i)
    {
        max = fabs(a[i * n + i]);
        max_i = i;
        for (int k = i + 1; k < n; ++k)
        {
            if (fabs(a[k * n + i]) > max)
            {
                max = fabs(a[k * n + i]);
                max_i = k;
            }
        }

        if (fabs(a[max_i * n + i]) < eps)
        {
            return false;
        }

        if (max_i != i)
        {
            for (j = 0; j < n; ++j)
            {
                SWAP(a[i * n + j], a[max_i * n + j]);
                SWAP(c[i * n + j], c[max_i * n + j]);
            }
        }

        d = a[i * n + i];
        for (j = 0; j < n; ++j)
        {
            a[i * n + j] /= d;
            c[i * n + j] /= d;
        }

        for (k = 0; k < n; ++k)
        {
            if (k != i)
            {
                d = a[k * n + i];
                for (j = 0; j < n; ++j)
                {
                    a[k * n + j] -= d * a[i * n + j];
                    c[k * n + j] -= d * c[i * n + j];
                }
            }
        }
    }
    return true;
}

void multy(double *a, double *b, int n, int p, int m, double *c)
{
    int i, j, k;
    double s, s00, s01, s02, s10, s11, s12, s20, s21, s22, a0, a1, a2, b0, b1, b2;

    int n3 = n - n % 3;
    int m3 = m - m % 3;

    for (i = 0; i < n3; i += 3)
    {
        for (j = 0; j < m3; j += 3)
        {
            s00 = s01 = s02 = s10 = s11 = s12 = s20 = s21 = s22 = 0;

            for (k = 0; k < p; k++)
            {
                a0 = a[i * p + k];
                a1 = a[(i + 1) * p + k];
                a2 = a[(i + 2) * p + k];

                b0 = b[k * m + j];
                b1 = b[k * m + j + 1];
                b2 = b[k * m + j + 2];

                s00 += a0 * b0;
                s01 += a0 * b1;
                s02 += a0 * b2;

                s10 += a1 * b0;
                s11 += a1 * b1;
                s12 += a1 * b2;

                s20 += a2 * b0;
                s21 += a2 * b1;
                s22 += a2 * b2;
            }
            c[i * m + j] = s00;
            c[i * m + j + 1] = s01;
            c[i * m + j + 2] = s02;

            c[(i + 1) * m + j] = s10;
            c[(i + 1) * m + j + 1] = s11;
            c[(i + 1) * m + j + 2] = s12;

            c[(i + 2) * m + j] = s20;
            c[(i + 2) * m + j + 1] = s21;
            c[(i + 2) * m + j + 2] = s22;
        }

        for (; j < m; j++)
        {
            s00 = s10 = s20 = 0;
            for (k = 0; k < p; k++)
            {
                s00 += a[i * p + k] * b[k * m + j];
                s10 += a[(i + 1) * p + k] * b[k * m + j];
                s20 += a[(i + 2) * p + k] * b[k * m + j];
            }
            c[i * m + j] = s00;
            c[(i + 1) * m + j] = s10;
            c[(i + 2) * m + j] = s20;
        }
    }
    for (; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            s = 0;
            for (k = 0; k < p; k++)
            {
                s += a[i * p + k] * b[k * m + j];
            }
            c[i * m + j] = s;
        }
    }
}

void multy_vector(double *g, double *c, int v, int h, double *d)
{
    int i, j;
    double s;
    for (i = 0; i < v; i++)
    {
        s = 0;
        for (j = 0; j < h; j++)
        {
            s += g[i * h + j] * c[j];
        }
        d[i] = s;
    }
}

void swap_rows(double *a, int n, int m, int p, int k, double *b, int s, int i0)
{
    int i, j;
    int s_loc = g2l_b(n, m, p, k, s);
    if (i0 != s_loc)
    {
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n - s * m; j++)
            {
                SWAP(a[n * (m * i0 + i) + m * s + j],
                     a[n * (m * s_loc + i) + m * s + j]);
            }
            SWAP(b[m * s_loc + i], b[m * i0 + i]);
        }
    }
}

void swap_columns(double *a, int n, int m, int p, int k, int s, int j0)
{
    int i, j;
    int rows = get_rows(n, m, p, k);
    // int s_loc = g2l_b(n, m, p, k, s);
    if (j0 != s)
    {
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < rows; i++)
            {
                SWAP(a[m * j0 + n * i + j],
                     a[m * s + n * i + j]);
            }
        }
    }
}

void sub(double *a, double *b, int n, int m, double *c)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            c[i * m + j] = a[i * m + j] - b[i * m + j];
        }
    }
}
int gaussian_method(double *a, double *b, double *x, int n, int m, int p, int k, double *buf, MPI_Comm com)
{
    int b_rows = get_block_rows(n, m, p, k), b_columns = (n + m - 1) / m;
    double *c = new double[m * m], *c_inv = new double[m * m], *d = new double[m * m], *e = new double[m * m], *f = new double[m * m],
           *buf_block_b = new double[m];
    int *permutation;
    int buf_block_b_size;
    int i, j;
    int v, h;
    int s;
    double norm_a = matrix_norm(a, n, m, p, k, com);

    double min_norm = 0, norm;
    int b_rows_m = (l2g_b(n, m, p, k, b_rows - 1) == n / m ? b_rows - 1 : b_rows);
    int b_columns_m = n / m;
    bool flag = true;
    int res_i, res_j;
    //   int start;

    permutation = new int[b_columns + 1];
    for (i = 0; i < b_columns + 1; i++)
    {
        permutation[i] = i;
    }
    for (s = 0; s < b_columns_m; s++)
    {
        int sk = s % p;
        int s_loc = g2l_b(n, m, p, k, s);

        // ---------------------------------- №1 ----------------------------------
        // границы для обхода по квадратным блокам
        // start = g2l_b(n, m, p, k, p * ((p - 1 + s - k) / p) + k);
        // printf("k = %d, start = %d\n", k, start);
        flag = true;
        i = 0;
        while (l2g_b(n, m, p, k, i) < s)
            i++;
        for (; i < b_rows_m; i++)
        {
            for (j = s; j < b_columns_m; j++)
            {
                get_block(a, n, m, p, k, i, j, c, v, h); // v = h = m
                // print_matrix_local(c, v, h, k);
                // if (k == 0)
                //     printf("\n");
                if (inverse(c, m, c_inv, norm_a))
                {
                    // print_matrix_local(c_inv, v, h, k);
                    norm = matrix_norm_local(c_inv, v);
                    // printf("k = %d, norm = %lf\n", k, norm);
                    if (flag)
                    {
                        flag = false;
                        min_norm = norm;
                        res_i = i;
                        res_j = j;
                    }
                    else if (min_norm > norm)
                    {
                        min_norm = norm;
                        res_i = i;
                        res_j = j;
                    }
                }
            }
        }
        // printf("Process: %d, min_norm = %lf, res_i = %d, res_j = %d\n", k, min_norm, res_i, res_j);
        //  находим глобальный минимум
        struct
        {
            double norm_inv; // обратная норма
            int res;         // i * n + j
        } local_pair, global_pair;
        if (flag) // все блоки вырождены
        {
            local_pair.norm_inv = -1;
            local_pair.res = -1;
        }
        else
        {
            local_pair.norm_inv = 1 / min_norm;
            local_pair.res = l2g_b(n, m, p, k, res_i) * b_columns + res_j;
        }
        MPI_Allreduce(&local_pair, &global_pair, 1, MPI_DOUBLE_INT, MPI_MAXLOC, com);
        if (global_pair.res == -1) // во всей матрице нет обратимых блоков
        {
            if (k == 0)
                printf("This method is not applicable with the given parameters\n");
            delete[] c;
            delete[] c_inv;
            delete[] d;
            delete[] e;
            delete[] f;
            delete[] buf_block_b;
            delete[] permutation;

            return 0;
        }
        // printf("k = %d, nrm_incv = %lf, res = %d\n", k, local_pair.norm_inv, local_pair.res);
        int res_i_glob = global_pair.res / b_columns, res_j_glob = global_pair.res % b_columns;
        res_i = g2l_b(n, m, p, k, res_i_glob);
        res_j = res_j_glob;

        // if (k == 0)
        // printf("Process: %d, global_min_norm = %lf, global_res_i = %d, global_res_j = %d\n", k, 1 / global_pair.norm_inv, res_i_glob, res_j_glob);
        // printf("l2_norm(a) = %lf\n", l2_norm_matrix(a, n, m, p, k, com));
        //  ---------------------------------- №2 ----------------------------------
        int owner = res_i_glob % p;
        if (owner == sk)
        {
            if (k == sk)
            {
                swap_rows(a, n, m, p, k, b, s, res_i);
            }
        }
        else
        {
            MPI_Status st;
            if (k == sk)
            {
                // memcpy(buf, a + 0 * n * m, n * m * sizeof(double));
                MPI_Sendrecv_replace(a + s_loc * n * m, n * m, MPI_DOUBLE, owner, 0, owner, 0, com, &st);
                MPI_Sendrecv_replace(b + s_loc * m, m, MPI_DOUBLE, owner, 0, owner, 0, com, &st);
            }
            else if (k == owner)
            {
                // memcpy(buf, a + res_i * n * m, n * m * sizeof(double));
                MPI_Sendrecv_replace(a + res_i * n * m, n * m, MPI_DOUBLE, sk, 0, sk, 0, com, &st);
                MPI_Sendrecv_replace(b + res_i * m, m, MPI_DOUBLE, sk, 0, sk, 0, com, &st);
            }
            // print_matrix_local(buf, m, n, k);
        }
        // printf("res_j = %d\n", res_j);
        swap_columns(a, n, m, p, k, s, res_j);

        std::swap(permutation[s], permutation[res_j]);

        // if (k == 0)
        //     printf("Matrix after selecting the main element %d %d, s = %d\n", res_i_glob, res_j_glob, s);
        // print_matrix(a, n, m, p, k, buf, MAX_PRINT, com);
        // if (k == 0)
        //     printf("\n");
        // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
        // if (k == 0)
        //     printf("\n");

        // printf("l2_norm(a) = %lf\n", l2_norm_matrix(a, n, m, p, k, com));

        //  ---------------------------------- №3 ----------------------------------
        //  домножение строки на обратную
        if (k == sk)
        {
            memcpy(buf, a + s_loc * n * m, n * m * sizeof(double));
        }
        MPI_Bcast(buf, n * m, MPI_DOUBLE, sk, com);

        int len = b_columns / p;
        int end = (k != p - 1 ? (k + 1) * len : b_columns);
        int w;
        get_block(buf, n, m, p, sk, 0, s, c, v, h); // v = h = m
        // if (k == 0)
        //     printf("C:\n");
        // print_matrix_local(c, m, m, k);
        inverse(c, m, c_inv, norm_a);
        // if (k == 0)
        //     printf("C_inv:\n");
        // print_matrix_local(c_inv, m, m, k);
        for (j = k * len; j < end; j++)
        {
            get_block_buf(buf, n, m, p, sk, s_loc, j, c, h, w);
            // print_matrix_local(c, h, w, k, p - 1);
            // printf("k = %d, v = %d,h = %d, w = %d\n", k, v, h, w);
            multy(c_inv, c, v, h, w, d);
            // print_matrix_local(d, v, w, k, p - 1);
            set_block(buf, n, m, 0, j, d, v, w);
        }
        if (k == sk)
        {
            get_block_vector(b, n, m, p, k, s_loc, c, h);
            // print_matrix_local(c, 1, h, k);
            multy_vector(c_inv, c, v, h, d);
            set_block_vector(b, m, s_loc, d, v);
            memcpy(buf_block_b, d, v * sizeof(double));
            buf_block_b_size = v;
        }
        MPI_Bcast(buf_block_b, v, MPI_DOUBLE, sk, com);
        MPI_Bcast(&buf_block_b_size, 1, MPI_INT, sk, com);
        // Обнуляем все то, что не трогал данный процесс
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < k * len * m; j++)
            {
                buf[i * n + j] = 0;
            }
        }
        for (i = 0; i < m; i++)
        {
            for (j = end * m; j < n; j++)
            {
                buf[i * n + j] = 0;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, buf, n * m, MPI_DOUBLE, MPI_SUM, com);
        if (k == sk)
        {
            memcpy(a + s_loc * n * m, buf, n * m * sizeof(double));
        }
        // printf("l2_norm(a) = %lf\n", l2_norm_matrix(a, n, m, p, k, com));

        // if (k == 0)
        //     printf("Matrix after row multiplication, s = %d\n", s);
        // print_matrix(a, n, m, p, k, buf, MAX_PRINT, com);
        // if (k == 0)
        //     printf("\n");
        // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
        // if (k == 0)
        //     printf("\n");
        // ---------------------------------- №4 ----------------------------------
        // преобразования строк
        int cv, ch, dv, dh, ev, eh;
        // int begin = g2l_b(n, m, p, k, p * ((p + s - k) / p) + k);
        //  if (k == 0)
        //      printf("k = %d, begin = %d\n", k, begin);
        //  print_matrix_local(buf_block_b, 1, buf_block_b_size, k);

        // устанавливаем на начало
        i = 0;
        while (l2g_b(n, m, p, k, i) <= s)
            i++;
        // printf(" k = %d, i = %d, s = %d\n", k, i, s);
        for (; i < b_rows; i++)
        {
            get_block(a, n, m, p, k, i, s, c, cv, ch);
            for (j = s; j < b_columns; j++)
            {
                get_block(a, n, m, p, k, i, j, d, dv, dh);
                get_block_buf(buf, n, m, p, sk, s_loc, j, e, ev, eh);
                multy(c, e, cv, ch, eh, f);
                sub(d, f, cv, eh, d);
                set_block(a, n, m, i, j, d, cv, eh);
            }
            get_block_vector(b, n, m, p, k, i, d, dv);
            // if (k == 0)
            //     printf("d:\n");
            // print_matrix_local(d, 1, dv, k);
            multy_vector(c, buf_block_b, cv, buf_block_b_size, f);
            // if (k == 0)
            //     printf("c:\n");
            // print_matrix_local(c, cv, ch, k);
            sub(d, f, 1, dv, d);
            set_block_vector(b, m, i, d, dv);
        }

        // if (k == 0)
        //     printf("Matrix after row transformation, s = %d\n", s);
        // print_matrix(a, n, m, p, k, buf, MAX_PRINT, com);
        // if (k == 0)
        //     printf("\n");
        // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
        // if (k == 0)
        //     printf("\n");
    }
    int sk = s % p;
    int s_loc = g2l_b(n, m, p, k, s);
    if (n % m != 0)
    {
        if (k == sk)
        {
            // printf("b_rows = %d, s_loc = %d\n", b_rows, s_loc);
            get_block(a, n, m, p, k, s_loc, s, c, v, h); // v = h = l
            // if (k == 0)
            //     printf("C:\n");
            //  print_matrix_local(c, v, h, k); // v = h
            if (!inverse(c, v, c_inv, norm_a))
            {
                if (k == 0)
                    printf("This method is not applicable with the given parameters\n");
                delete[] c;
                delete[] c_inv;
                delete[] d;
                delete[] e;
                delete[] f;
                delete[] buf_block_b;
                delete[] permutation;

                return 0;
            }
            // можно убрать
            for (i = 0; i < v; i++)
            {
                for (j = 0; j < h; j++)
                {
                    e[i * h + j] = (i == j ? 1 : 0);
                }
            }
            set_block(a, n, m, s_loc, s, e, v, h);

            get_block_vector(b, n, m, p, k, s_loc, c, h);
            // print_matrix_local(c, 1, h, k);
            multy_vector(c_inv, c, v, h, d);
            set_block_vector(b, m, s_loc, d, v);
        }
    }
    // if (k == 0)
    //     printf("Matrix after main step, s = %d\n", s);
    // print_matrix(a, n, m, p, k, buf, MAX_PRINT, com);
    // if (k == 0)
    //     printf("\n");
    // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
    // if (k == 0)
    //     printf("\n");

    int cv, ch, dv;
    // if (n % m != 0)
    // {
    //     MPI_Bcast(d, m * m, MPI_DOUBLE, sk, com);
    //     dv = v;
    // }
    // else
    // {
    //     sk = (s - 1) % p;
    //     MPI_Bcast(d, m * m, MPI_DOUBLE, sk, com);
    // }

    for (s = b_columns - 1; s >= 0; s--)
    {
        int owner = s % p;
        int s_loc = g2l_b(n, m, p, k, s);
        if (k == owner)
        {
            get_block_vector(b, n, m, p, k, s_loc, d, dv);
        }
        MPI_Bcast(&dv, 1, MPI_DOUBLE, owner, com);
        MPI_Bcast(d, dv, MPI_DOUBLE, owner, com);
        // print_matrix_local(d, 1, dv, k);
        for (i = 0; i < b_rows && l2g_b(n, m, p, k, i) < s; i++)
        {
            get_block(a, n, m, p, k, i, s, c, cv, ch);

            // if (k == 0)
            //     printf("a[i,s]:\n");
            // print_matrix_local(c, cv, ch, k);

            multy_vector(c, d, cv, ch, e);
            get_block_vector(b, n, m, p, k, i, c, cv);

            // if (k == 0)
            //     printf("b[i]:\n");
            // print_matrix_local(c, 1, cv, k);

            sub(c, e, 1, cv, c);
            set_block_vector(b, m, i, c, cv);
        }
        // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
    }
    // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
    // for (i = 0; i < b_rows; i++)
    // {
    //     int i_glob = l2g_b(n, m, p, k, i);
    //     // меняем i_glob c permutation[i_glob]
    //     int pk = permutation[i_glob] % p;
    //     int perm_loc = g2l_b(n, m, p, k, permutation[i_glob]);
    //     printf("k = %d, pk = %d, i = %d, perm_loc = %d \n", k, pk, i, perm_loc);
    //     if (k == pk)
    //     {
    //         // swap_rows(a, n, m, p, k, b, i, perm_loc);
    //         memcpy(x, b + perm_loc * m, m * sizeof(double));
    //         memcpy(b + perm_loc * m, b + i * m, m * sizeof(double));
    //         memcpy(b + i * m, x, m * sizeof(double));
    //     }
    //     // else
    //     // {
    //     //     MPI_Status st;
    //     //     MPI_Sendrecv_replace(b + i * m, m, MPI_DOUBLE, k, 0, k, 0, com, &st);
    //     //     if (k == pk)
    //     //     {
    //     //         MPI_Sendrecv_replace(b + perm_loc * m, m, MPI_DOUBLE, pk, 0, pk, 0, com, &st);
    //     //     }
    //     //     // print_matrix_local(buf, m, n, k);
    //     // }
    // }

    for (i = 0; i < b_columns; i++)
    {

        int owner_pi = permutation[i] % p, owner_i = i % p;
        int perm_loc = g2l_b(n, m, p, k, permutation[i]);
        int i_loc = g2l_b(n, m, p, k, i);
        // printf("perm_loc = %d, i_loc = %d\n", perm_loc, i_loc);
        if (owner_pi == owner_i)
        {
            if (k == owner_pi)
            {
                memcpy(x + perm_loc * m, b + i_loc * m, m * sizeof(double));
                //  memcpy(b + perm_loc * m, b + i_loc * m, m * sizeof(double));
                // memcpy(b + i_loc * m, x, m * sizeof(double));
            }
        }
        else
        {
            MPI_Status st;
            if (k == owner_pi)
            {
                MPI_Recv(x + perm_loc * m, m, MPI_DOUBLE, owner_i, 0, com, &st);
            }
            else if (k == owner_i)
            {
                MPI_Send(b + i_loc * m, m, MPI_DOUBLE, owner_pi, 0, com);
            }
        }
    }

    // if (k == 0)
    // {
    //     for (i = 0; i < b_columns; i++)
    //     {
    //         printf(" %d", permutation[i]);
    //     }
    //     printf("\n");
    // }

    // if (k == 0)
    //     printf("Matrix after all, s = %d\n", s);
    // print_matrix(a, n, m, p, k, buf, MAX_PRINT, com);
    // if (k == 0)
    //     printf("\n");
    // print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
    // if (k == 0)
    //     printf("\n");

    // memcpy(x, b, b_rows * m * sizeof(double));
    delete[] c;
    delete[] c_inv;
    delete[] d;
    delete[] e;
    delete[] f;
    delete[] buf_block_b;
    delete[] permutation;
    return 1;
}