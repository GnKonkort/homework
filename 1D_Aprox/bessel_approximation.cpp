#include "bessel_approximation.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
bessel_approximation::bessel_approximation(/* args */)
{
}

bessel_approximation::~bessel_approximation()
{
}
double bessel_approximation::separated_differences(double y1, double y2, double x1, double x2){
    return (y2 - y1) / (x2 - x1);
}
double bessel_approximation::d(int i,int n,double* x, double* y, double* y_derivatives){
    if((i == 0 || i == n-1)){
        return y_derivatives[i];
    } else {
        return (y[i+1] - y[i - 1])/(x[i+1] - x[i-1]);
    }
}
double* bessel_approximation::build_coeffecients(int n, double *x, double *y, double *y_derivatives){
    clock_t tStart = clock();
    double* result = new double[4*(n)];
    // a_{1,i} = f[x_i]
    // a_{2,i} = d_i
    // a_{3,i} = (f(x_i,x_{i+1}) - d_i) / (x_{i+1} - x_i)
    // a_{4,i} = (d_i + d_{i+1} - 2*f(x_i,x_{i+1})) / (x_{i+1} - x_i)

    for(int i = 0; i < n; i++){
//        result[i*4] = y[i];
//        result[i*4 + 1] = d(i,n,x,y,y_derivatives);
//        result[i*4 + 2] = (separated_differences(y[i],y[i+1],x[i],x[i+1]) - d(i,n,x,y,y_derivatives)) / (x[i+1] - x[i]);
//        result[i*4 + 3] = (d(i,n,x,y,y_derivatives) + d(i+1,n,x,y,y_derivatives) - 2 * separated_differences(y[i],y[i+1],x[i],x[i+1])) / ((x[i+1] - x[i]) * (x[i+1] - x[i]));
        result[i*4] = y[i];
        result[i*4 + 1] = d(i,n,x,y,y_derivatives);
        result[i*4 + 2] = (3*separated_differences(y[i],y[i+1],x[i],x[i+1]) - 2 * d(i,n,x,y,y_derivatives) - d(i+1,n,x,y,y_derivatives)) / (x[i+1] - x[i]);
        result[i*4 + 3] = (d(i,n,x,y,y_derivatives) + d(i+1,n,x,y,y_derivatives) - 2 * separated_differences(y[i],y[i+1],x[i],x[i+1])) / ((x[i+1] - x[i]) * (x[i+1] - x[i]));
    }
    //printf("[DEBUG] Table built for %d * 4 = %d points,it took %.10fs seconds\n",n , 4*n, (double)(clock() - tStart) / CLOCKS_PER_SEC);
    fflush_unlocked(stdout);
    return result;
}



double bessel_approximation::approx_function(double x,int n,double *x_massive ,double* table){
    int counter = 0;
    double h = x_massive[1] - x_massive[0];
    counter = (int)(fabs(x - x_massive[0]) / h);
    if(counter < 0){
            counter = 0;
        }
        if(counter >= n){
            counter = n;
        }
    //printf("[DEBUG] %lf belongs to [%lf,%lf]\n",x,x_massive[counter],x_massive[counter+1]);
    return table[counter*4] + table[counter*4 + 1]*(x - x_massive[counter]) + table[counter*4 + 2]*(x - x_massive[counter])*(x - x_massive[counter]) + table[counter*4 + 3]*(x - x_massive[counter])*(x - x_massive[counter])*(x - x_massive[counter]);
}
