#ifndef __THREAD_F_H__
#define __THREAD_F_H__
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
void *thread_function(void *g);
struct args{
    pthread_mutex_t* mutex;
    pthread_barrier_t* barrier;
    pthread_cond_t *condvar;

    int n,m,r,s,l;
    int k,p;
    int error;
    double *A,*B;
    double norm;
    double *reverse_block_buffer;
    double reverse_block_norm;
    int reverse_block_pos;
    double elapsed_thread;
    double wall_clock;
    args* a;
};

int solve(double *A, double *B, int n, int m, int k, int p, pthread_mutex_t* mutex,pthread_barrier_t* barrier);

#endif