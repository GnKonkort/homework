#ifndef BESSEL_APPROXIMATION_H
#define BESSEL_APPROXIMATION_H

class bessel_approximation
{
private:

public:
    bessel_approximation();
    ~bessel_approximation();
    double* build_coeffecients(int n, double *x, double *y, double *y_derivatives);
    double d(int n,int i,double* x, double* y, double* y_derivatives);
    double separated_differences(double y1, double y2, double x1, double x2);
    double approx_function(double x,int n,double *x_massive ,double* table);
};




#endif //BESSEL_APPROXIMATION
