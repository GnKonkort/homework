#include "newtonaproximator.h"
#include <math.h>
#include <iostream>
double NewtonAproximator::derivative(double x, double(*func)(double)){
    double h = 1.e-6;
    return (func(x + h) - func(x - h))/(2*h);
}

void NewtonAproximator::build_table(int n, double (*func)(double),double a, double b, int mistake_amount){
    this->n = n;
    this->a = a;
    this->b = b;
    this->mistake_amount = mistake_amount;
    this->func = func;
    table.reset(new double[2*n]);
    points.reset(new double[2*n]);
    double dx = (b - a) / (n-1);
    fmax = fabs(func(a));
    for(int i = 0; i < n; i++){
        if(fmax < fabs(func(a + i*dx))) fmax = fabs(func(a + i*dx));
    }
    for(int i = 0; i < n; i++){
        points[2*i] = a + i*dx;
        points[2*i + 1] = a + i*dx;
        table[2*i] = local_func(a + i*dx);
        table[2*i + 1] = local_func(a + i*dx);
    }
    //First iteration
    for(int i = 2*n - 1; i > 0; i--){
        if(fabs(points[i] - points[i - 1]) < 1.e-6){
            table[i] = derivative(points[i],func);
        } else {
            table[i] = (table[i] - table[i-1]) / (points[i] - points[i-1]);
        }
    }
    //Every concurrent iteration
    for(int i = 1; i < 2*n; i++){
        for(int j = 2*n - 1; j > i; j--){
            table[j] = (table[j] - table[j - 1]) / (points[j] - points[j - i - 1]);
        }
    }
}

double NewtonAproximator::aprox(double x){
    double result = 0.0;
    double buffer = 0.0;
    for(int i = 0; i < 2*n; i++){
        buffer = 1.0;
        for(int j = 0; j < i; j++){
            buffer *= (x - points[j]);
            //printf("[DEBUG]Buffer = %.3e\n",buffer);
        }
        buffer *= table[i];
        result += buffer;
    }
    return result;
}

double NewtonAproximator::local_func(double x){
    if(fabs(x - (b + a) / (2)) < 1.e-6){
        return func(x) + fmax * 0.1 * mistake_amount;
    }
    return func(x);
}

NewtonAproximator::NewtonAproximator(){
}
NewtonAproximator::~NewtonAproximator(){
    points.release();
    table.release();
}
