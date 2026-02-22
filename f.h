#ifndef F_H
#define F_H
#include "mpi.h"

int g2l_b(int n, int m, int p, int k, int i_glob);
int l2g_b(int n, int m, int p, int k, int i_loc);
int l2g(int n, int m, int p, int k, int i_loc);
int get_block_rows(int n, int m, int p, int k);
int get_rows(int n, int m, int p, int k);
int get_max_block_rows(int n, int m, int p, int k);
double f(int s, int n, int i, int j);

void init_matrix(double *a, int n, int m, int p, int k, int s);
int read_array(FILE *fp, double *a, int n);
int read_matrix(double *a, int n, int m, int p, int k, const char *name, double *buf, MPI_Comm com);
void print_matrix(double *a, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm com);
int print_array(double *a, int n, int m, int printed_rows, int max_print);
void print_vector(double *b, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm com);
int print_vector_part(double *v, int rows, int printed, int max_print);
void matrix_mult_vector(double *a, double *b, double *c, int n, int m, int p, int k, MPI_Comm com);
double l1_norm(double *x, int n, int m, int p, int k, MPI_Comm com);

void init_b(double *a, double *b, int n, int m, int p, int k);
void init_x_exact(double *x, int n, int m, int p, int k);

double get_r1(double *a, double *x, double *b, double *buf, int n, int m, int p, int k, MPI_Comm com);
double get_r2(double *x, double *x_exact, int n, int m, int p, int k, MPI_Comm com);

void get_block(double *a, int n, int m, int p, int k, int i, int j, double *c, int &v, int &h);
void set_block(double *a, int n, int m, int i, int j, double *c, int v, int h);
void print_matrix_local(double *c, int v, int h, int k);
void gaussian_method(double *a, double *b, double *x, int n, int m, int p, int k, MPI_Comm com);
#endif