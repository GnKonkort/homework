#ifndef THREAD_KERNEL_H
#define THREAD_KERNEL_H
#include <pthread.h>
#include "aproximator.h"
class thread_kernel
{
private:
    pthread_t* threads;
    ThreadArgument* arguments;
    int nx,ny;
public:
    Databus* common;
    static void* thread_function(void* a);
    void initialize_threads();
    void terminate_threads();
    double f_local(double x, double y);
    void assign_new_task(int nx, int ny,double x0,double x1, double y0, double y1, double (*f)(double, double), int error_amount, double eps);
    thread_kernel(int p);
};

#endif // THREAD_KERNEL_H
