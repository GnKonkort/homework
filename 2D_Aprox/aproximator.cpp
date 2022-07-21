#include "aproximator.h"
#include <cstring>
#include <math.h>

double aproximator::f(double x, double y)
{
    double a = (this->common->x1 + this->common->x0) / this->common->nx * (this->common->nx / 2);
    double b = (this->common->y1 + this->common->y0) / this->common->ny * (this->common->ny / 2);
    if(fabs(x - a) < __DBL_EPSILON__ && fabs(y - b) < __DBL_EPSILON__){
        return this->common->f(x,y) + 0.1 * this->common->error_amount * this->common->f_max;
    }
    return this->common->f(x,y);
}

void aproximator::reduce_sum(int p, int *a, int n)
{
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t cin = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t cout = PTHREAD_COND_INITIALIZER;

    static int tin = 0;
    static int tout = 0;
    static int *p_a = 0;
    int i = 0;

    if(p <= 1)
        return;

    pthread_mutex_lock(&mutex);
    if(!p_a) p_a = a;
    else
        for(i=0; i<n; i++)
            p_a[i] += a[i];

    tin++;
    if(tin >= p){
        tout = 0;
        pthread_cond_broadcast(&cin);
    }
    else
        while(tin < p)
            pthread_cond_wait(&cin, &mutex);

    if(p_a != a)
        for(i=0; i<n; i++)
            a[i] = p_a[i];

    tout++;
    if(tout >= p){
        tin = 0;
        p_a = 0;
        pthread_cond_broadcast(&cout);
    }
    else
        while(tout < p)
            pthread_cond_wait(&cout, &mutex);

    pthread_mutex_unlock(&mutex);
}

int aproximator::get_len_msr(int nx, int ny)
{
    return 6*(nx-1)*(ny-1) + 4*((nx-1)*2 + (ny-1)*2) + 3*2 + 2*2;
}

void aproximator::ij2l(int nx, int ny, int i, int j, int &l)
{
    (void) ny;
    l = i + j*(nx+1);
}

void aproximator::l2ij(int nx, int ny, int &i, int &j, int l)
{
    (void) ny;
    j = l/(nx+1);
    i = l - j*(nx+1);
}

void aproximator::precond(double *A, int n, double *y, double *x, int p, int k)
{
    //printf("[DEBUG] x = %lf y = %lf A = %lf\n",x[0],y[0],A[0]);
    ////fflush(stdout);
    int i1, i2, i;

    if(n%p == 0){
        i1 = k*(n/p);
        i2 = (k+1)*(n/p);
    }
    else{
        if(k < n%p){
            i1 = k*(n/p) + k;
            i2 = (k+1)*(n/p) + k+1;
        }
        else{
            i1 = n%p + k*(n/p);
            i2 = n%p + (k+1)*(n/p);
        }
    }

    for(i=i1; i<i2; i++){
        //printf("[%d]\n",i);
        ////fflush(stdout);
        x[i] = y[i]/A[i];
    }
    reduce_sum(p, 0, 1);
}

double aproximator::scalar_product(int n, double *x, double *y, double *buf, int p, int k)
{
    int i1, i2;
    double s = 0;
    int i;

    if(n%p == 0){
        i1 = k*(n/p);
        i2 = (k+1)*(n/p);
    }
    else{
        if(k < n%p){
            i1 = k*(n/p) + k;
            i2 = (k+1)*(n/p) + k+1;
        }
        else{
            i1 = n%p + k*(n/p);
            i2 = n%p + (k+1)*(n/p);
        }
    }

    for(i=i1; i<i2; i++){
        s += x[i]*y[i];
    }
    buf[k] = s;
    reduce_sum(p, 0, 1);

    s = 0;
    for(i=0; i<p; i++) s += buf[i];
    reduce_sum(p, 0, 1);
    return s;
}

void aproximator::mult_scal_plus_vector(double *x, double *y, int n, double t, int p, int k)
{
    int i1, i2, i;

    if(n%p == 0){
        i1 = k*(n/p);
        i2 = (k+1)*(n/p);
    }
    else{
        if(k < n%p){
            i1 = k*(n/p) + k;
            i2 = (k+1)*(n/p) + k+1;
        }
        else{
            i1 = n%p + k*(n/p);
            i2 = n%p + (k+1)*(n/p);
        }
    }

    for(i=i1; i<i2; i++){
        x[i] -= t*y[i];
    }
    reduce_sum(p, 0, 1);
}

void aproximator::allocate_arrays(int nx, int ny, int p, int* &I, double* &A, double* &B, double* &U, double* &X, double* &R, double* &V, double* &buf)
{
    int tmp = (nx + 1)*(ny + 1);

    if(I != nullptr) delete[] I;
    if(A != nullptr) delete[] A;
    if(B != nullptr) delete[] B;
    if(U != nullptr) delete[] U;
    if(X != nullptr) delete[] X;
    if(R != nullptr) delete[] R;
    if(V != nullptr) delete[] V;
    if(buf != nullptr) delete[] buf;



    I = new int[tmp + 1 + get_len_msr(nx, ny)];
    A = new double[tmp + 1 + get_len_msr(nx, ny)];
    B = new double[tmp];
    U = new double[tmp];
    X = new double[tmp];
    R = new double[tmp];
    V = new double[tmp];
    buf = new double[p];

    memset(B, 0., tmp*sizeof(double));
    memset(U, 0., tmp*sizeof(double));
    memset(X, 0., tmp*sizeof(double));
    memset(R, 0., tmp*sizeof(double));
    memset(V, 0., tmp*sizeof(double));
    memset(A, 0., (tmp + 1 + get_len_msr(nx, ny))*sizeof(double));
    memset(I, 0, (tmp + 1 + get_len_msr(nx, ny))*sizeof(int));
    memset(buf, 0, p*sizeof(double));
}

int aproximator::allocate_msr_matrix(int nx, int ny, double *A, int *I)
{
    int N = (nx+1)*(ny+1);
    int nz = get_len_msr(nx, ny);
    int len = N+1+nz;

    //A = new double[len];
    if(!A) return -1;
    //I = new int[len];
    if(!I){
        delete []A;
        return -1;
    }

    I[0] = N+1;
    int k, l=0, s=0;
    for(k=1; k<=N; k++){
        l = get_off_diag_num(nx, ny, k-1);
        I[k] = I[k-1] + l;
        //printf("I[%d] = %d\n", k, I[k]);
        s += l;
    }

    if(s!=nz){
        printf ("[FATAL] Cannot allocate MSR\n");
        return -1;
    }
    if(I[N]!=len){
        printf ("[FATAL] Cannot allocate MSR\n");
        return -1;
    }

    return 0;
}

void aproximator::build_msr_matrix(int nx, int ny, double *A, int *I, int p, int k)
{
    int k1, k2, i;
    int s=0, sum=0;;
    int N = (nx+1)*(ny+1);

    if(N%p == 0){
        k1 = k*(N/p);
        k2 = (k+1)*(N/p);
    }
    else{
        if(k < N%p){
            k1 = k*(N/p) + k;
            k2 = (k+1)*(N/p) + k+1;
        }
        else{
            k1 = N%p + k*(N/p);
            k2 = N%p + (k+1)*(N/p);
        }
    }

    for(i=k1; i<k2; i++){
        s = get_off_diag_index(nx, ny, i, I+I[i], A+i, A+I[i]);
        sum += s;
    }

    A[N] = 0.;
    reduce_sum(p, &sum, 1);

    if((sum + N + 1) != I[N]){
        printf ("Error InitA\n");
    }
}

void aproximator::allocate_msr_vec(double *b, int nx, int ny, int p, int k/*, double (*f)(double, double)*/, double x0, double y0, double x1, double y1)
{
    int i1 = 0, i2 = 0, i = 0;
    int N = (nx+1)*(ny+1);
    if(N%p == 0){
        i1 = k*(N/p);
        i2 = (k+1)*(N/p);
    }
    else{
        if(k < N%p){
            i1 = k*(N/p) + k;
            i2 = (k+1)*(N/p) + k+1;
        }
        else{
            i1 = N%p + k*(N/p);
            i2= N%p + (k+1)*(N/p);
        }
    }

    for(i=i1; i<i2; i++){
        b[i] = integrate_msr(nx, ny, i/*, f*/, x0, y0, x1, y1);
        //printf("[DEBUG] i1 = %d i2 = %d b[%d] = %lf\n",i1,i2,i,b[i]);
        ////fflush(stdout);
    }
    reduce_sum(p, 0, 1);
}

double aproximator::integrate_msr(int nx, int ny, int l/*, double (*f)(double, double)*/, double x2, double y2, double x1, double y1)
{
    int i = 0, j = 0;

    l2ij(nx, ny, i, j, l);

    double hx = (x1 - x2)/nx;
    double hy = (y2 - y1)/ny;
    double x0 = x2+i*hx, y0 = y1+j*hy;
    double l1 = 0, l2 = 0, l3 = 0, l4 = 0;
    //printf("[AAAAA] %lf\n",f(10,10));
    if(i>0 && i<nx && j>0 && j<ny){
        l1 = f(x0, y0);
        l2 = f(x0 + hx/2, y0) + f(x0 + hx/2, y0 + hy/2) + f(x0, y0 + hy/2) + f(x0 - hx/2, y0) + f(x0 - hx/2, y0 - hy/2) + f(x0, y0 - hy/2);
        l3 = f(x0 + hx, y0 + hy/2) +  f(x0 + hx/2, y0 + hy) + f(x0 - hx/2, y0 + hy/2) + f(x0 - hx, y0 - hy/2) + f(x0 - hx/2, y0 - hy) + f(x0 + hx/2, y0 - hy/2);
        l4 = f(x0 + hx, y0) + f(x0 + hx, y0 + hy) + f(x0, y0 + hy) + f(x0 - hx, y0) + f(x0 - hx, y0 - hy) + f(x0, y0 - hy);

        return (18*l1 + 10*l2 + 2*l3 + l4)/96;
    }

    if(i==0 && j>0 && j<ny){
        l1 = f(x0, y0);
        l2 = 2*f(x0 + hx/2, y0) + 2*f(x0 + hx/2, y0 + hy/2) + f(x0, y0 + hy/2) + f(x0, y0 - hy/2);
        l3 = f(x0 + hx, y0 + hy/2) + f(x0 + hx/2, y0 + hy) + f(x0 + hx/2, y0 - hy/2);
        l4 = 2*f(x0 + hx, y0) + 2*f(x0 + hx, y0 + hy) + f(x0, y0 + hy) + f(x0, y0 - hy);

        return (18*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i==nx && j>0 && j<ny){
        l1 = f(x0, y0);
        l2 = f(x0, y0 + hy/2) + 2*f(x0 - hx/2, y0) + 2*f(x0 - hx/2, y0 - hy/2) + f(x0, y0 - hy/2);
        l3 = f(x0 - hx/2, y0 + hy/2) + f(x0 - hx, y0 - hy/2) + f(x0 - hx/2, y0 - hy);
        l4 = f(x0, y0 + hy) + 2*f(x0 - hx, y0) + 2*f(x0 - hx, y0 - hy) + f(x0, y0 - hy);

        return (18*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i>0 && i<nx && j==0){
        l1 = f(x0, y0);
        l2 = f(x0 + hx/2, y0) + 2*f(x0 + hx/2, y0 + hy/2) + 2*f(x0, y0 + hy/2) + f(x0 - hx/2, y0);
        l3 = f(x0 + hx, y0 + hy/2) + f(x0 + hx/2, y0 + hy) + f(x0 - hx / 2, y0 + hy/2);
        l4 = f(x0 + hx, y0) + 2*f(x0 + hx, y0 + hy) + 2*f(x0, y0 + hy) + f(x0 - hx, y0);

        return (18*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i>0 && i<nx && j==ny){
        l1 = f(x0, y0);
        l2 = f(x0 + hx/2, y0) + f(x0 - hx/2, y0) + 2*f(x0 - hx/2, y0 - hy/2) + 2*f(x0, y0 - hy/2);
        l3 = f(x0 - hx, y0 - hy/2) + f(x0 - hx/2, y0 - hy) + f(x0 + hx/2, y0 - hy/2);
        l4 = f(x0 - hx, y0) + 2*f(x0 - hx, y0 - hy) + 2*f(x0, y0 - hy) + f(x0 + hx, y0);

        return (18*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i==0 && j==0){
        l1 = f(x0, y0);
        l2 = f(x0 + hx/2, y0) + 2*f(x0 + hx/2, y0 + hy/2) + f(x0, y0 + hy/2);
        l3 = f(x0 + hx, y0 + hy/2) + f(x0 + hx/2, y0 + hy);
        l4 = f(x0 + hx, y0) + 2*f(x0 + hx, y0 + hy) + f(x0, y0 + hy);

        return (12*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i==nx && j==ny){
        l1 = f(x0, y0);
        l2 = f(x0 - hx/2, y0) + 2*f(x0 - hx/2, y0 - hy/2) + f(x0, y0 - hy/2);
        l3 = f(x0 - hx, y0 - hy/2) + f(x0 - hx/2, y0 - hy);
        l4 = f(x0 - hx, y0) + 2*f(x0 - hx, y0 - hy) + f(x0, y0 - hy);

        return (12*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i==0 && j==ny){
        l1 = f(x0, y0);
        l2 = f(x0 + hx/2, y0) + f(x0, y0 - hy/2);
        l3 = f(x0 + hx/2, y0 - hy/2);
        l4 = f(x0 + hx, y0) + f(x0, y0 - hy);

        return (6*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    if(i==nx && j==0){
        l1 = f(x0, y0);
        l2 = f(x0, y0 + hy/2) + f(x0 - hx/2, y0);
        l3 = f(x0 - hx/2, y0 + hy/2);
        l4 = f(x0, y0 + hy) + f(x0 - hx, y0);

        return (6*l1 + 10*l2 + 4*l3 + l4)/192;
    }

    return -1;
}

void aproximator::mat_vec_mul_msr(int n, double *A, int *I, double *x, double *y, int p, int k)
{
    int i, j, l;
    double s;
    int J;

    int i1, i2;
    if(n%p == 0){
        i1 = k*(n/p);
        i2 = (k+1)*(n/p);
    }
    else{
        if(k < n%p){
            i1 = k*(n/p) + k;
            i2 = (k+1)*(n/p) + k+1;
        }
        else{
            i1 = n%p + k*(n/p);
            i2 = n%p + (k+1)*(n/p);
        }
    }

    for(i=i1; i<i2; i++){
        s = A[i]*x[i];     //y_i = a_ii*x_i
        l = I[i+1] - I[i]; //число не равных 0 элементво i-й строки
        J = I[i];        //начало строки i не равной 0, внедиагональный элемент
        for(j=0; j<l; j++){
            s += A[J+j]*x[I[J+j]];
            y[i] = s;
        }
    }
    reduce_sum(p, 0, 1);
}

int aproximator::get_off_diag_index(int nx, int ny, int l, int *J, double *a_diag, double *a)
{
    int i, j;
    l2ij(nx, ny, i, j, l);

    if(i>0 && i<nx && j>0 && j<ny){
        ij2l(nx, ny, i+1,   j, J[0]);
        ij2l(nx, ny, i,   j-1, J[1]);
        ij2l(nx, ny, i-1, j-1, J[2]);
        ij2l(nx, ny, i-1,   j, J[3]);
        ij2l(nx, ny, i,   j+1, J[4]);
        ij2l(nx, ny, i+1, j+1, J[5]);

        *a_diag = 1./2;
        for(int ii=0; ii<6; ii++){
            a[ii] = 1./12;
        }

        return 6;
    }

    if(i==0 && j>0 && j<ny){ //левая
        ij2l(nx, ny, i+1,   j, J[0]);
        ij2l(nx, ny, i,   j-1, J[1]);
        ij2l(nx, ny, i,   j+1, J[2]);
        ij2l(nx, ny, i+1, j+1, J[3]);

        *a_diag = 1./4;
        a[0] = 1./12;
        a[1] = 1./24;
        a[2] = 1./24;
        a[3] = 1./12;

        return 4;
    }

    if(i==nx && j>0 && j<ny){ //правая
        ij2l(nx, ny, i,   j-1, J[0]);
        ij2l(nx, ny, i-1, j-1, J[1]);
        ij2l(nx, ny, i-1,   j, J[2]);
        ij2l(nx, ny, i,   j+1, J[3]);

        *a_diag = 1./4;
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./12;
        a[3] = 1./24;

        return 4;
    }

    if(i>0 && i<nx && j==0){ //низ
        ij2l(nx, ny, i+1,   j, J[0]);
        ij2l(nx, ny, i-1,   j, J[1]);
        ij2l(nx, ny, i,   j+1, J[2]);
        ij2l(nx, ny, i+1, j+1, J[3]);

        *a_diag = 1./4;
        a[0] = 1./24;
        a[1] = 1./24;
        a[2] = 1./12;
        a[3] = 1./12;

        return 4;
    }

    if(i>0 && i<nx && j==ny){ //верх
        ij2l(nx, ny, i+1,   j, J[0]);
        ij2l(nx, ny, i,   j-1, J[1]);
        ij2l(nx, ny, i-1, j-1, J[2]);
        ij2l(nx, ny, i-1,   j, J[3]);

        *a_diag = 1./4;
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./12;
        a[3] = 1./24;

        return 4;
    }

    if(i==0 && j==0){
        ij2l(nx, ny, i+1,   j, J[0]);
        ij2l(nx, ny, i,   j+1, J[1]);
        ij2l(nx, ny, i+1, j+1, J[2]);

        *a_diag = 1./6;
        a[0] = 1./24;
        a[1] = 1./24;
        a[2] = 1./12;

        return 3;
    }

    if(i==nx && j==ny){
        ij2l(nx, ny, i,   j-1, J[0]);
        ij2l(nx, ny, i-1, j-1, J[1]);
        ij2l(nx, ny, i-1,   j, J[2]);

        *a_diag = 1./6;
        a[0]=1./24;
        a[1]=1./12;
        a[2]=1./24;

        return 3;
    }

    if(i==0 && j==ny){
        ij2l(nx, ny, i+1,   j, J[0]);
        ij2l(nx, ny, i,   j-1, J[1]);

        *a_diag = 1./12;
        a[0] = 1./24;
        a[1] = 1./24;

        return 2;
    }

    if(i==nx && j==0){
        ij2l(nx, ny, i-1,   j, J[0]);
        ij2l(nx, ny, i,   j+1, J[1]);

        *a_diag = 1./12;
        a[0] = 1./24;
        a[1] = 1./24;

        return 2;
    }

    return -1;
}

int aproximator::get_off_diag_num(int nx, int ny, int l)
{
    int i, j;
    l2ij(nx, ny, i, j, l);

    if(i>0 && i<nx && j>0 && j<ny) return 6;

    if(i==0 && j>0 && j<ny) return 4;  //левая

    if(i==nx && j>0 && j<ny) return 4;  //правая

    if(i>0 && i<nx && j==0) return 4;  //низ

    if(i>0 && i<nx && j==ny) return 4;  //верх

    if(i==0 && j==0) return 3;

    if(i==nx && j==ny) return 3;

    if(i==0 && j==ny) return 2;

    if(i==nx && j==0) return 2;

    return -1;
}

int aproximator::solver(int n, double *A, int *I, double *b, double *x, double *r, double *u, double *v, double eps, int maxit, int p, int k, double *buf, int *iterations)
{
    int step = 10;
    int ret=0, it=0;
    for(it=0; it<maxit; it+=step){
        //printf("[DEBUG SOLVER] it = %d\n",it);
        ret += solver_iteration(n, A, I, b, x, r, u, v, eps, step, p, k, buf);
        //printf("[DEBUG] Current ret = %d\n",ret);
        ////fflush(stdout);
        reduce_sum(p, 0, 1);
        if(ret >= 0) break;
    }

    iterations[0] = ret;
    //printf("[DEBUG] Amogus\n");
    //////fflush(stdout);
    reduce_sum(p, 0, 1);
    //for(int iii=0; iii<10; iii++) printf("x[%d] = %lf\n", iii, x[iii]);
    if(k==0) printf("[APROXIMATOR] Finished with number of iterations: %d\n", iterations[0]);
    if(it >= maxit) return -1;
    ////fflush(stdout);
    return 0;
}

int aproximator::solver_iteration(int n, double *A, int *I, double *b, double *x, double *r, double* u, double *v, double eps, int maxit, int p, int k, double *buf)
{
    int it;
    double prec, b_norm2;
    double C1, C2, tau;

    b_norm2 = scalar_product(n, b, b, buf, p, k);
    //printf("[DEBUG] b_norm2 = %.3e eps = %.3e\n",b_norm2,eps);
    prec = eps*eps*b_norm2;
    // r = Ax
    mat_vec_mul_msr(n, A, I, x, r, p, k);
    // r -= b
    mult_scal_plus_vector(r, b, n, 1., p, k);
    //printf("[DEBUG] prec = %.3e\n",prec);
    for(it=0; it<maxit; it++){

        precond(A, n, r, v, p, k);

        mat_vec_mul_msr(n, A, I, v, u, p, k);
        ////fflush(stdout);
        C1 = scalar_product(n, u, r, buf, p, k);
        C2 = scalar_product(n, u, u, buf, p, k);
        ////fflush(stdout);

        if(C1 < prec || C2 < prec)
            break;
        ////fflush(stdout);
        tau = C1/C2;
        ////fflush(stdout);
        mult_scal_plus_vector(x, v, n, tau, p, k);
        mult_scal_plus_vector(r, u, n, tau, p, k);
    }

    reduce_sum(p, 0, 1);
    if(it >= maxit) return -1;
    //printf("[DEBUG] Solver iteration sucessful, returning with code %d\n",it);
    ////fflush(stdout);
    return it;
}

void aproximator::normalize_msr(double *A, double *b, int nx, int ny, int p, int k, double x0, double y0, double x1, double y1, double LEN, int N)
{
    (void)p;
    (void)k;
    int i = 0;
    //printf("N = %d\n", N);
    A[N] = 0.;
    //printf("%lf, %lf\n %lf, %lf\n\n", x0, y0, x1, y1);
    //printf("A[0] = %lf\n", A[0]);
    for(i=0; i<N+LEN+1; i++)
        A[i] = A[i]*(x1 - x0)*(y0 - y1)/(nx*ny);
    for(i=0; i<N; i++)
        b[i] = b[i]*(x1 - x0)*(y0 - y1)/(nx*ny);
}

aproximator::aproximator(Databus *common)
{
    this->common = common;
}



double Databus::f_local(double x, double y)
{
    double a = (this->x1 + this->x0) / nx * (nx / 2);
    double b = (this->y1 + this->y0) / nx * (nx / 2);
    if(fabs(x - a) < __DBL_EPSILON__ && fabs(y - b) < __DBL_EPSILON__){
        return this->f_max * 0.1 * this->error_amount + this->f(x,y);

    }
    return f(x,y);
}
