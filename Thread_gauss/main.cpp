#include <stdio.h>
#include <stdlib.h>
#include "libs/matrix.h"
#include "libs/thread_function.h"
#include <limits>
#include <fenv.h>
int main(int argc, char** argv){
    //a.out n m p r s [filename]
    feenableexcept(FE_UNDERFLOW || FE_OVERFLOW || FE_INVALID || FE_DIVBYZERO);
    int n,m,p,r,s;
    char filename[256];
    if(argc > 7 || argc < 6){
        printf("Usage: %s n m p r s [filename]\n",argv[0]);
        return -1;
    }
    if( !sscanf(argv[1],"%d",&n)||
        !sscanf(argv[2],"%d",&m)||
        !sscanf(argv[3],"%d",&p)||
        !sscanf(argv[4],"%d",&r)||
        !sscanf(argv[5],"%d",&s)){
            printf("Usage: %s n m p r s [filename]\n",argv[0]);
            return -1;
        }
    if(s == 0){
        if(argc != 7){
            printf("Usage: %s n m p r s [filename]\n",argv[0]);
            return -1;
        }
        if(!sscanf(argv[6],"%s",filename)){
            printf("Usage: %s n m p r s [filename]\n",argv[0]);
            return -1;
        }
    }
    if(n < 0 ||
       m > n ||
       p <= 0 ||
       s > 4||
       s < 0){
        printf("Usage: %s n m p r s [filename]\n",argv[0]);
        return -1;
    }
    std::unique_ptr<double[]> A{new double[n*n]},B{new double[n*n]},C{new double[n*n]},backup_matrix{new double[n*n]},reverse_block_buffer{new double[m*m]};
    if(s != 0){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                A[i*n + j] = backup_matrix[i*n + j] = formula(n,i,j,s);
            }
        }
    } else {
        if(read_matrix(n,0,A.get(),B.get(),fopen(filename,"r")) == -1){
            printf("Usage: %s n m p r s [filename]\n",argv[0]);
            return -1;
        }
        
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            B[i*n + j] = (i == j);
        }
        for(int i = 0; i < n*n; i++){
            backup_matrix[i] = A[i];
        }
    }
    args* arguments = new args[p];
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier,0,p);
    for(int i = 0; i < p; i++){
        arguments[i].k = i;
        arguments[i].p = p;
        arguments[i].n = n;
        arguments[i].m = m;
        arguments[i].r = r;
        arguments[i].s = s;
        arguments[i].A = A.get();
        arguments[i].B = B.get();
        arguments[i].reverse_block_buffer = reverse_block_buffer.get();
        arguments[i].barrier = &barrier;
        arguments[i].mutex = &mutex;
        arguments[i].a = arguments;
        arguments[i].norm = norma(A.get(),n,n);
        arguments[i].reverse_block_norm = std::numeric_limits<double>::max();
        arguments[i].reverse_block_pos = -1;
    }
    for(int i = 0; i < p; i++){
        printf("[DEBUG] reverse_block_buffer pointer %d = %p\n",i,arguments[i].reverse_block_buffer);
    }
    pthread_t* threads = new pthread_t[p];
    printf("[DEBUG] A:\n");
    print_matrix(A.get(),n,n,r);
    printf("[DEBUG] B:\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    print_matrix(B.get(),n,n,r);
    for (int k = 1; k < p; k++)
    if (pthread_create(threads + k, 0, &thread_function, arguments + k))
      printf("Cannot create thread %d\n", k);
    thread_function(arguments + 0);

    for (int k = 1; k < p; k++) pthread_join(threads[k], 0);
    pthread_barrier_destroy(&barrier);
    for(int i = 0; i < p; i++){
        if (arguments[i].error == -1){
            printf("Matrix is irreversable\n");
            return -1;
        }
    }
    printf("[DEBUG] A:\n");
    print_matrix(A.get(),n,n,r);
    printf("[DEBUG] B:\n");
    print_matrix(B.get(),n,n,r);
    double elapsed = 0.0;
    for(int i = 0; i < p; i++){
        if(arguments[i].wall_clock > elapsed){
            elapsed = arguments[i].wall_clock;
        }
    }
    multiplication(B.get(),backup_matrix.get(),C.get(),n,n,n);
    for(int i = 0; i < n; i++){
        C[i*n + i] -= 1.0;
    }
    printf("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], norma(C.get(),n,n), elapsed, s, n, m);

    return 0;
}