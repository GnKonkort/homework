#ifndef NEWTONAPROXIMATOR_H
#define NEWTONAPROXIMATOR_H


#include <memory>
class NewtonAproximator{
    private:
    int n;
    double a,b;
    int mistake_amount;
    double fmax;
    double (*func)(double);
    std::unique_ptr<double[]> table;
    std::unique_ptr<double[]> points;

    public:
    NewtonAproximator();
    ~NewtonAproximator();
    void build_table(int n, double (*func)(double),double a, double b,int mistake_amount);
    double local_func(double x);
    double aprox(double x);
    double derivative(double x, double(*func)(double));
};
#endif // NEWTONAPROXIMATOR_H
