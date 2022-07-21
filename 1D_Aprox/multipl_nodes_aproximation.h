#ifndef MULTIPL_NODES_APROXIMATION_H
#define MULTIPL_NODES_APROXIMATION_H


class Multipl_nodes_aproximation
{
public:
    Multipl_nodes_aproximation();
    //Функиця построения таблицы для схемы Горвера с кратными узлами
    double* build_table(int n,double *x, double *y,double *y_derivatives);
    double approximation_function(int n,double x, double* x_massive,double *table);
};

#endif // MULTIPL_NODES_APROXIMATION_H

