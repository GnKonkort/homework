#ifndef __MPI_MATRIX_H__
#define __MPI_MATRIX_H__
#include <math.h>
#include <mpich/mpi.h>
#include <cstring>
#include <memory>
#include <limits>
int get_max_rows(int n, int m, int p);
double formula_1(int i, int j, int n);
double formula_2(int i, int j, int n);
double formula_3(int i, int j, int n);
double formula_4(int i, int j, int n);
double formula_E(int i, int j, int n);
double formula_0(int i, int j, int n);
int get_rows(int n, int m, int p, int k);
int get_rows_block(int n, int m, int p, int k);
void init_matrix(double *a, int n, int m, int p, int k,double(*f)(int i, int j, int n));
int read_matrix(double *a, int n, int m, int p, int k, const char *name, double *buf /* n*m */);
void print_matrix(const double *a, int n, int m, int p, int k, double *buf, int max_print);
int read_array(FILE* fp,double* a, int len);
double parallel_matrix_norm(double* a,int k, int rows, int n);
int solve(double *A, double *B, int n, int m, int p, int k, double* buffer_a, double* buffer_b, double norm);
void mpi_matrix_multiplication(double *a, double *b,double *c, int n, int m, int p, int k);
#endif