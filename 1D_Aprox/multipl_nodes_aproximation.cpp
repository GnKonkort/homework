#include "multipl_nodes_aproximation.h"
#include <math.h>
#include <stdio.h>
Multipl_nodes_aproximation::Multipl_nodes_aproximation()
{

}
double* Multipl_nodes_aproximation::build_table(int n,double *x, double *y, double *y_derivatives){
    //Выделение памяти для таблицы
    double* table = new double[2*n];
    //Первая итерация заполнения матрицы
    for(int i = 0; i < 2*n; i++){
        table[i] = y[i];
    }
    //Вторая итерация заполнения (отдельно выделена, поскольку тут нужно знать производные)
    for(int i = 2*n - 1; i >= 1; i--){
        //printf("[DEBUG] Comparing %lf <-> %lf\n",table[i],table[i - 1]);
        if(fabs(table[i] - table[i - 1]) < __DBL_EPSILON__){
            table[i] = y_derivatives[i];
        } else {
            table[i] = (table[i] - table[i - 1]) / (x[i] - x[i - 1]);
        }
    }
    //Все последующие итерации
    for(int i = 2; i < 2*n; i++){
        for(int j = 2*n - 1; j >= i; j--){
            //printf("[debug] table[j] - table[j-1] / x[j] - x[j-i] = %lf - %lf / %lf - %lf=",table[j],table[j-1],x[j],x[j-i]);
            fflush(stdout);
            table[j] = (table[j] - table[j - 1]) / (x[j] - x[j - i]);
            //printf("%lf\n",table[j]);
        }
    }
    return table;
}
double Multipl_nodes_aproximation::approximation_function(int n,double x, double* x_massive,double *table){
    double s = table[2*n - 1];
        for (int i = 2; i <= 2*n; ++i)
        {
            s *= (x - x_massive[2*n - i]);
            s += table[2*n-i];
        }
        return s;
}
