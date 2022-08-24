#include "besselaproximator.h"
#include <math.h>
double BesselAproximator::derivative(double x, double(*func)(double)){
    double h = 1.e-6;
    return (func(x + h) - func(x - h))/(2*h);
}

double BesselAproximator::d(double x, double (*func)(double), double dx, double a, double b){
    if(fabs(x - a) < 1.e-6 || fabs(x - b) < 1.e-6){
        return derivative(x,func);
    }
    return (local_func(x + dx) - local_func(x - dx)) / (2*dx);
}


void BesselAproximator::build_table(int n, double a, double b, double (*func)(double), int mistake_amount){
    this->n = n;
    this->a = a;
    this->b = b;
    this->func = func;
    this->dx = (b - a) / (n - 1);
    this->mistake_amount = mistake_amount;
    this->fmax = fabs(func(a));
    for(int i = 0; i < n; i++){
        if(fabs(func(a + i*dx)) > this->fmax) this->fmax = this->func(a + i * dx);
    }
    c1.reset(new double[n]);
    c2.reset(new double[n]);
    c3.reset(new double[n]);
    c4.reset(new double[n]);
    double x = 0;
    for(int i = 0; i < n; i++){
        x = a + dx * i;
        c1[i] = local_func(x);
        c2[i] = d(x,func,dx,a,b);
        c3[i] = ((local_func(x + dx) - local_func(x))/(dx) - d(x,func,dx,a,b))/
                (dx);
        c4[i] = (d(x,func,dx,a,b) + d(x + dx, func,dx,a,b) - 2 * (local_func(x + dx) - local_func(x))/(dx)) /
                (dx*dx);
    }
}

double BesselAproximator::aprox(double x){
    int i = static_cast<int>(std::floor((x - a) / dx));
    //printf("[DEBUG] %lf is in [%lf,%lf] = %d\n",x,a + dx*i, a + dx*(i+1),i);
    return c1[i] + (x - (a + dx*i))*c2[i] + (x - (a + dx*i))*(x - (a + dx*i))*c3[i] + (x - (a + dx*i))*(x - (a + dx*i))*(x - (a + dx*i))*c4[i];
}

double BesselAproximator::local_func(double x){
    if(fabs(x - (b + a) / 2) < 1.e-6){
        return func(x) + 0.1 * fmax * mistake_amount;
    }
    return func(x);
}

BesselAproximator::BesselAproximator(){

}
BesselAproximator::~BesselAproximator(){
    c1.release();
    c2.release();
    c3.release();
    c4.release();
}
