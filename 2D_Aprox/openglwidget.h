#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H
#include <QOpenGLWidget>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QString>
#include <QFont>
#include "thread_kernel.h"
#include <QTimer>
#include "grid.h"

enum GraphMode : int{
    REAL_FUNC,
    APROX_FUNC,
    RESIDUAL
};
enum FuncMode : int{
    CONST,
    X,
    Y,
    SUM,
    SQUARED_ROOT,
    SQUARED,
    EXPONENT,
    NON_POLYNOMIAL
};


class openglwidget : public QOpenGLWidget
{
    Q_OBJECT
private:
    int nx,ny;
    int p;
    double x0,y0,x1,y1;
    double angle_psi, angle_fi;
    double scale;
    GraphMode graph_mode;
    FuncMode func_mode;
    double (*f)(double,double);
    double f_local(double x, double y);
    QTimer* thread_overwatch;
    double* X_local;
    bool changes_blocked;
    int error_amount;
    double f_min,f_max;
    double residual;
    double eps;
public:
    Grid* painting_grid;
    thread_kernel* thread_controller;
    QSize sizeHint();
    void initializeGL();
    void paintGL();
    openglwidget(QWidget* parent, int nx, int ny, double x0, double x1, double y0, double y1, int p, int start_function,double eps);
    void fill_grid();
public slots:
    void rotate_cw();
    void rotate_ccw();
    void double_scale();
    void sub_scale();
    void double_nx();
    void sub_nx();
    void double_ny();
    void double_nx_ny();
    void sub_nx_ny();
    void sub_ny();
    void change_func();
    void change_graph();
    void check_threads();
    void increment_error();
    void decrement_error();
};

#endif // OPENGLWIDGET_H
