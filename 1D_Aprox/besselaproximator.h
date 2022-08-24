#ifndef BESSELAPROXIMATOR_H
#define BESSELAPROXIMATOR_H


#include <memory>
class BesselAproximator{
    private:
    int n;
    double a,b;
    double dx;
    double (*func)(double);
    double fmax;
    int mistake_amount;
    std::unique_ptr<double[]> c1,c2,c3,c4;
    public:
    void build_table(int n, double a, double b, double (*func)(double), int mistake_amount);
    double local_func(double x);
    double aprox(double x);
    double d(double x, double (*func)(double), double dx, double a, double b);
    double derivative(double x, double(*func)(double));
    BesselAproximator();
    ~BesselAproximator();
};

#endif // BESSELAPROXIMATOR_H
