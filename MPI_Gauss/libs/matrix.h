#ifndef __MATRIX_H__
#define __MATRIX_H__


#include <memory>
int inverse_matrix(double *a, double *x, int n, double norm);
void print_matrix(double* A, int n, int m,int r);
void get_block(double* A, int n, int m, int i, int j, double* block);
void put_block(double* A, int n, int m, int i, int j, double* block);
double formula(int n, int x, int y, int mode);
double norma(double *a, int p, int q);
void multiplication(double *a, double *b, double *c, int n, int m,int q);
void ones(double *a, int n);
void substraction(double *a, double *b, int n, int m);
void zero(double *a, int n, int m);
void get_block_raw(double* A, int n, int m, int i, int j, double* block, int p, int q);
void put_block_raw(double* A, int n, int m, int i, int j, double* block, int p, int q);
#endif