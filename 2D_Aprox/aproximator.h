#ifndef APROXIMATOR_H
#define APROXIMATOR_H
#include <stdio.h>
#include <pthread.h>

enum ThreadStatus{
    FREE,
    BUSY
};

enum TaskStatus{
    NO_TASK,
    CHANGED_NXNY,
    CHANGED_FUNC,
    TASK_COMPLETE
};

//This class allows public access to all vactors/matrices
class Databus{
public:
    int *I = nullptr;
    double *A = nullptr;
    double *B = nullptr;
    double *X = nullptr;
    double *U = nullptr;
    double *V = nullptr;
    double *R = nullptr;
    double *buf = nullptr;
    double eps;
    double f_max;
    int nx,ny, LEN, iterations;
    TaskStatus current_task_status;
    int error_amount;
    double x0,y0,x1,y1;
    double (*f)(double x, double y);
    double f_local(double x, double y);
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mutex_calculation = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

    Databus() = default;
    ~Databus() = default;
};

class ThreadArgument{
public:
    int iterations;
    ThreadStatus thread_status;
    Databus* common;
    int p,k;
};


//This class provides all information required for threads to perform aproximation
class aproximator
{
public:
    Databus* common;
    double (*f_real)(double x, double y);
    double f(double x, double y);
    //thread synchronizer
    void reduce_sum(int p, int *a, int n);

    //Coordinates transformation (i,j) -> l and l ->(i,j)
    int get_len_msr(int nx, int ny);
    void ij2l(int nx, int ny, int i, int j, int &l);
    void l2ij(int nx, int ny, int &i, int &j, int l);

    //Jacobi's preconditioner for residual minimization
    void precond(double *A, int n, double *y, double *x, int p, int k);

    //Usual matrix operations(but parallel)
    double scalar_product(int n, double *x, double *y, double* buffer, int p, int k);
    void mult_scal_plus_vector(double *x, double*y, int n, double t, int p, int k);

    //memory allocations
    void allocate_arrays(int nx, int ny, int p, int* &I, double* &A, double* &B, double* &U, double* &X, double* &R, double* &V, double* &buf);
    //preparations for MSR matrix and MSR vector(in previous version these functions were responsible for memory allocation, now they only check and fill)
    int allocate_msr_matrix(int nx, int ny, double *A, int *I);
    void build_msr_matrix(int nx, int ny, double *A, int *I, int p, int k); // separate function to fill matrix with values
    void allocate_msr_vec(double *b, int nx, int ny, int p, int k/*, double (*f)(double, double)*/, double x0, double y0, double x1, double y1);
    double integrate_msr(int nx, int ny, int l/*, double (*f)(double, double)*/, double x2, double y2, double x1, double y1); // get msr basis function integrations with current function
    void mat_vec_mul_msr(int n, double *A, int *I, double *x, double *y, int p, int k);
    //function for filling MSR matrix
    int get_off_diag_index(int nx, int ny, int l, int *J, double *a_diag, double *a);
    //function for filling index massive I
    int get_off_diag_num(int nx, int ny, int l);

    //solver function and solver iteration
    int solver(int n, double *A, int *I, double *b, double *x, double *r, double *u, double *v, double eps, int maxit, int p, int k, double *buf, int *iterations);
    int solver_iteration(int n, double *A, int *I, double *b, double *x, double *r, double* u, double *v, double eps, int maxit, int p, int k, double *buf);

    //MSR matrix normalization
    void normalize_msr(double *A, double *b, int nx, int ny, int p, int k, double x0, double y0, double x1, double y1, double LEN, int N);
    aproximator(Databus* common);
    ~aproximator() = default;
};

#endif // APROXIMATOR_H
