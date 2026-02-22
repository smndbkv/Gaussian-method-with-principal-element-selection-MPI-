#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <sys/sysinfo.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include "io_status.h"

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
    double sum_loc = 0, sum_glob = 0;
    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        sum_loc += fabs(x[i_loc] - x_exact[i_loc]);
    }
    MPI_Allreduce(&sum_loc, &sum_glob, 1, MPI_DOUBLE, MPI_SUM, com);
    return sum_glob;
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

void print_matrix_local(double *c, int v, int h, int k)
{
    for (int i = 0; i < v; i++)
    {
        for (int j = 0; j < h; j++)
        {
            if (k == 0)
            {
                printf(" %10.3e", c[i * h + j]);
            }
        }
        if (k == 0)
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
                SWAP(a[n * (m * i0 + i) + m * s_loc + j],
                     a[n * (m * s_loc + i) + m * s_loc + j]);
            }
            SWAP(b[m * s_loc + i], b[m * i0 + i]);
        }
    }
}

void swap_columns(double *a, int n, int m, int p, int k, int s, int j0)
{
    int i, j;
    int rows = get_rows(n, m, p, k);
    int s_loc = g2l_b(n, m, p, k, s);
    if (j0 != s)
    {
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < rows; i++)
            {
                SWAP(a[m * j0 + n * i + j],
                     a[m * s_loc + n * i + j]);
            }
        }
    }
}

void gaussian_method(double *a, double *b, double *x, int n, int m, int p, int k, double *buf, MPI_Comm com)
{
    int b_rows = get_block_rows(n, m, p, k), b_columns = (n + m - 1) / m;
    double *c = new double[m * m], *c_inv = new double[m * m];
    int i, j;
    int v, h;
    double norm_a = matrix_norm(a, n, m, p, k, com);

    double min_norm = 0, min_norm_glob, norm;
    int b_rows_m = (l2g_b(n, m, p, k, b_rows - 1) == n / m ? b_rows - 1 : b_rows);
    int b_columns_m = n / m;
    bool flag = true;
    int res_i, res_j;
    int start;
    for (int s = 0; s < 2; s++)
    {
        // ---------------------------------- №1 ----------------------------------
        // ищем квадратный блок с минимальной нормой обратной матрицы
        // границы для обхода по квадратным блокам
        start = g2l_b(n, m, p, k, p * ((p - 1 + s - k) / p) + k);
        // printf("k = %d, start = %d\n", k, start);
        flag = true;
        for (i = start; i < b_rows_m; i++)
        {
            for (j = s; j < b_columns_m; j++)
            {
                get_block(a, n, m, p, k, i, j, c, v, h); // v = h
                // print_matrix_local(c, v, h, k);
                // if (k == 0)
                //     printf("\n");
                if (inverse(c, m, c_inv, norm_a))
                {
                    print_matrix_local(c_inv, v, h, k);
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
            return;
        }
        int res_i_glob = global_pair.res / b_columns, res_j_glob = global_pair.res % b_columns;
        res_i = g2l_b(n, m, p, k, res_i_glob);
        res_j = res_j_glob;

        // printf("Process: %d, global_min_norm = %lf, global_res_i = %d, global_res_j = %d\n", k, 1 / global_pair.norm_inv, res_i, res_j);

        // ---------- №2 ----------
        // переставляем строки и столбцы
        int owner = res_i_glob % p;
        int sk = s % p;
        int s_loc = g2l_b(n, m, p, k, s);
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

        if (k == 0)
            printf("\n");
        print_matrix(a, n, m, p, k, buf, MAX_PRINT, com);
        if (k == 0)
            printf("\n");
        print_vector(b, n, m, p, k, buf, MAX_PRINT, com);
    }
    delete[] c;
    delete[] c_inv;
}