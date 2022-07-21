
#include <QPainter>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "window.h"
#include "multipl_nodes_aproximation.h"
#include "bessel_approximation.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10

static
double f_0(double x){
    return 1;
}

static
double f_1(double x){
    return x;
}
static
double f_2(double x){
    return x*x;
}
static
double f_3(double x){
    return x*x*x;
}
static
double f_4(double x){
    return x*x*x*x;
}
static
double f_5(double x){
    return exp(x);
}
static
double f_6(double x){
    return 1 / (25 * x * x + 1);
}
static
double derivative(double (*f)(double), double x){
    return (f(x + 1e-6) - f(x - 1e-6)) / (2 * 1e-6);
}
Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;
  f_name = "f (x) = 1";
  f = f_0;
  func_id = 0;
  approx_id = 0;
  mistake_scale = 0;
//  change_func ();
//  change_approx();
  update();
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
    return -1;

  if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
      || n <= 0)
    return -2;

  return 0;
}

/// change current function for drawing
void Window::change_func ()
{
  func_id = (func_id + 1) % 7;

  switch (func_id)
    {
      case 0:
        f_name = "f (x) = 1";
        f = f_0;
        break;
      case 1:
        f_name = "f (x) = x";
        f = f_1;
        break;
      case 2:
        f_name = "f (x) = x * x";
        f = f_2;
        break;
      case 3:
        f_name = "f (x) = x * x * x";
        f = f_3;
        break;
      case 4:
        f_name = "f (x) = x * x * x * x";
        f = f_4;
        break;
      case 5:
        f_name = "f (x) = exp(x)";
        f = f_5;
        break;
      case 6:
        f_name = "f (x) = 1 / (25*x*x + 1)";
        f = f_6;
        break;
    }
  update ();
}

void Window::double_n(){
    this->n *= 2;
    update();
}
void Window::sub_n(){
    if((this->n / 2) != 1){
        this->n /= 2;
    }
    update();
}

void Window::add_mistake()
{
    mistake_scale++;
    update();
}

void Window::sub_mistake()
{
    mistake_scale--;
    update();
}
void Window::change_approx(){
    approx_id = (approx_id + 1) % 4;
    switch (approx_id) {
        case 0:
            approx_name = "Multiple nodes approximation";
            break;
        case 1:
            approx_name = "Bessel polynom approximation";
            break;
        case 2:
            approx_name = "Multiple nodes and polynom approximation";
            break;
        case 3:
            approx_name = "Residual";
            break;
    }
    update();
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{
  double residual = 0.0;
  QSize widget_size = this->size();
  double aprox_painter_delta = (b - a) / (widget_size.rwidth());
  QPainter painter (this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  double *x, *y, *y_derivatives, *table;
  double delta_y, delta_x = (b - a) / (n-1);
  QPen pen_black(Qt::black, 0, Qt::SolidLine); 
  QPen pen_red(Qt::red, 0, Qt::SolidLine); 
  QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
  QPen pen_green(Qt::green, 0, Qt::SolidLine);
  clock_t tStart;

  Multipl_nodes_aproximation mul_node_approx;
  bessel_approximation bess_aprox;
  painter.setPen (pen_black);

  // calculate min and max for current function
  max_y = min_y = 0;
  for (x1 = a; x1 - b < 1.e-6; x1 += delta_x)
    {
      y1 = f (x1);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;
  printf("[DEBUG] Current mistake amount: %lf\n",this->mistake_scale);
  fflush_unlocked(stdout);
  // save current Coordinate System
  painter.save ();

  // make Coordinate Transformations
  painter.translate (0.5 * width (), 0.5 * height ());
  painter.scale (width () / (b - a), -height () / (max_y - min_y));
  painter.translate (-0.5 * (a + b), -0.5 * (min_y + max_y));
  painter.setPen (pen_red);
  painter.drawLine (a, 0, b, 0);
  painter.drawLine (0, max_y, 0, min_y);
  // draw approximated line for graph
    x1 = a;
    y1 = f (x1);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
      {
        y2 = f (x2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

        x1 = x2, y1 = y2;
      }
    x2 = b;
    y2 = f (x2);
    painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

  // draw approximated line for aproximation function
  switch (approx_id) {
    case 0:

      x = new double[2*n];

      y = new double[2*n];
      y_derivatives = new double[2*n];

      for(int i = 0; i < n; i++){
          x[2*i] = a + delta_x * i;
          x[2*i + 1] = a + delta_x * i;
          y[2*i] = f(a + delta_x * i);
          y[2*i + 1] = f(a + delta_x * i);
          if(i == n / 2){
              y[2*i] += max_y * 0.1 * mistake_scale;
              y[2*i + 1] += max_y * 0.1 * mistake_scale;
          }
          y_derivatives[2*i] = derivative(f,a + delta_x * i);
          y_derivatives[2*i + 1] = derivative(f,a + delta_x * i);
      }
//      printf("[DEBUG] Adding %lf * %lf * %d = %lf\n",max_y,0.1,mistake_scale,max_y * 0.1 * mistake_scale);
//      y[n / 2] += max_y * 0.1 * mistake_scale;
//      y[n / 2 + 1] += max_y * 0.1 * mistake_scale;
      painter.setPen (pen_blue);
      table = mul_node_approx.build_table(n,x,y,y_derivatives);

      x1 = a;
      y1 = mul_node_approx.approximation_function(n,x1,x,table);
      if(fabs(y1 - f(x1)) > residual){
          residual = fabs(y1 - f(x1));
      }
      //printf("[DEBUG] Point [%lf,%lf]\n",x1,y1);
      for (x2 = x1 + aprox_painter_delta; x2 - b < 1.e-6; x2 += aprox_painter_delta)
        {
          y2 = mul_node_approx.approximation_function(n,x2,x,table);
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          if(fabs(y2 - f(x2)) > residual){
              residual = fabs(y2 - f(x2));
          }
          //printf("[DEBUG] Point [%lf,%lf]\n",x2,y2);
          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = mul_node_approx.approximation_function(n,x2,x,table);
      if(fabs(y2 - f(x2)) > residual){
          residual = fabs(y2 - f(x2));
      }
      printf("[RESIDUAL]:%3e\n",residual);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
//      painter.setPen ("blue");
//      painter.drawText (0, 40, approx_name);
//      free(x);
//      free(y);
//      free(y_derivatives);
//      free(table);
      break;
    case 1:
      tStart = clock();
      x =  new double[n];
      y = new double[n];
      y_derivatives = new double[n];
      for(int i = 0; i < n; i++){
          x[i] = a + delta_x * i;
          y[i] = f(a + delta_x * i);
          y_derivatives[i] = derivative(f,a + delta_x * i);
      }
      //y_derivatives[0] = derivative(f,a);
      //y_derivatives[n-1] = derivative(f,a + delta_x * (n-1));
      painter.setPen (pen_blue);
      y[n / 2] += max_y * mistake_scale;
      table = bess_aprox.build_coeffecients(n,x,y,y_derivatives);
      //printf("[DEBUG] Table built for %d * 4 = %d points,it took %.10fs seconds\n",n , 4*n, (double)(clock() - tStart) / CLOCKS_PER_SEC);
      //tStart = clock();
      x1 = a;
      y1 = bess_aprox.approx_function(x1,n,x,table);
      if(fabs(y1 - f(x1)) > residual){
          residual = fabs(y1 - f(x1));
      }
      //printf("[DEBUG] Point [%lf,%lf]\n",x1,y1);
      for (x2 = x1 + aprox_painter_delta; x2 - b < 1.e-6; x2 += aprox_painter_delta)
        {
          y2 = bess_aprox.approx_function(x2,n,x,table);
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          if(fabs(y2 - f(x2)) > residual){
              residual = fabs(y2 - f(x2));
          }
          //printf("[DEBUG] Point [%lf,%lf]\n",x2,y2);
          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = bess_aprox.approx_function(x2,n,x,table);
      if(fabs(y2 - f(x2)) > residual){
          residual = fabs(y2 - f(x2));
      }
      printf("[RESIDUAL]:%3e\n",residual);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
//      painter.setPen ("blue");
//      painter.drawText (0, 40, approx_name);
//      free(x);
//      free(y);
//      free(y_derivatives);
//      free(table);
      break;
     case 2:
      x = new double[2*n];

      y = new double[2*n];
      y_derivatives = new double[2*n];
      painter.setPen(pen_blue);
      for(int i = 0; i < n; i++){
          x[2*i] = a + delta_x * i;
          x[2*i + 1] = a + delta_x * i;
          y[2*i] = f(a + delta_x * i);
          y[2*i + 1] = f(a + delta_x * i);
          y_derivatives[2*i] = derivative(f,a + delta_x * i);
          y_derivatives[2*i + 1] = derivative(f,a + delta_x * i);
      }
      painter.setPen (pen_blue);
      y[n / 2] += max_y * mistake_scale;
      table = mul_node_approx.build_table(n,x,y,y_derivatives);

      x1 = a;
      y1 = mul_node_approx.approximation_function(n,x1,x,table);
      if(fabs(y1 - f(x1)) > residual){
          residual = fabs(y1 - f(x1));
      }
      //printf("[DEBUG] Point [%lf,%lf]\n",x1,y1);
      for (x2 = x1 + aprox_painter_delta; x2 - b < 1.e-6; x2 += aprox_painter_delta)
        {
          y2 = mul_node_approx.approximation_function(n,x2,x,table);
          if(fabs(y2 - f(x2)) > residual){
              residual = fabs(y2 - f(x2));
          }
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          //printf("[DEBUG] Point [%lf,%lf]\n",x2,y2);
          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = mul_node_approx.approximation_function(n,x2,x,table);
      if(fabs(y2 - f(x2)) > residual){
          residual = fabs(y2 - f(x2));
      }
      printf("[RESIDUAL]:%3e\n",residual);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
//      painter.setPen ("blue");
//      painter.drawText (0, 40, approx_name);
      free(x);
      free(y);
      free(y_derivatives);
      free(table);
      painter.setPen(pen_green);
      x =  new double[n];
      y = new double[n];
      y_derivatives = new double[n];
      for(int i = 0; i < n; i++){
          x[i] = a + delta_x * i;
          y[i] = f(a + delta_x * i);
          y_derivatives[i] = derivative(f,a + delta_x * i);
      }
      //y_derivatives[0] = derivative(f,a);
      //y_derivatives[n-1] = derivative(f,a + delta_x * (n-1));
      painter.setPen (pen_green);
      y[n / 2] += max_y * mistake_scale;
      table = bess_aprox.build_coeffecients(n,x,y,y_derivatives);
      //printf("[DEBUG] Table built for %d * 4 = %d points,it took %.10fs seconds\n",n , 4*n, (double)(clock() - tStart) / CLOCKS_PER_SEC);
      //tStart = clock();
      x1 = a;
      y1 = bess_aprox.approx_function(x1,n,x,table);
      if(fabs(y1 - f(x1)) > residual){
          residual = fabs(y1 - f(x1));
      }
      //printf("[DEBUG] Point [%lf,%lf]\n",x1,y1);
      for (x2 = x1 + aprox_painter_delta; x2 - b < 1.e-6; x2 += aprox_painter_delta)
        {
          y2 = bess_aprox.approx_function(x2,n,x,table);
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          if(fabs(y2 - f(x2)) > residual){
              residual = fabs(y2 - f(x2));
          }
          //printf("[DEBUG] Point [%lf,%lf]\n",x2,y2);
          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = bess_aprox.approx_function(x2,n,x,table);
      printf("[RESIDUAL]:%3e\n",residual);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
//      painter.setPen ("blue");
//      painter.drawText (0, 40, approx_name);
//      free(x);
//      free(y);
//      free(y_derivatives);
//      free(table);
      break;
     case 3:
      x = new double[2*n];

      y = new double[2*n];
      y_derivatives = new double[2*n];

      for(int i = 0; i < n; i++){
          x[2*i] = a + delta_x * i;
          x[2*i + 1] = a + delta_x * i;
          y[2*i] = f(a + delta_x * i);
          y[2*i + 1] = f(a + delta_x * i);
          y_derivatives[2*i] = derivative(f,a + delta_x * i);
          y_derivatives[2*i + 1] = derivative(f,a + delta_x * i);
      }
      painter.setPen (pen_blue);
      table = mul_node_approx.build_table(n,x,y,y_derivatives);

      x1 = a;
      y1 = fabs(f(x1) - mul_node_approx.approximation_function(n,x1,x,table));
      //printf("[DEBUG] Point [%lf,%lf]\n",x1,y1);
      for (x2 = x1 + aprox_painter_delta; x2 - b < 1.e-6; x2 += aprox_painter_delta)
        {
          y2 = fabs(f(x2) - mul_node_approx.approximation_function(n,x2,x,table));
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          //printf("[DEBUG] Point [%lf,%lf]\n",x2,y2);
          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = fabs(f(x2) - mul_node_approx.approximation_function(n,x2,x,table));
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
//      painter.setPen ("blue");
//      painter.drawText (0, 40, approx_name);
//      free(x);
//      free(y);
//      free(y_derivatives);
//      free(table);
      x =  new double[n];
      y = new double[n];
      y_derivatives = new double[n];
      for(int i = 0; i < n; i++){
          x[i] = a + delta_x * i;
          y[i] = f(a + delta_x * i);
          y_derivatives[i] = derivative(f,a + delta_x * i);
      }
      painter.setPen (pen_green);
      table = bess_aprox.build_coeffecients(n,x,y,y_derivatives);
      x1 = a;
      y1 = fabs(f(x1) - bess_aprox.approx_function(x1,n,x,table));
      //printf("[DEBUG] Point [%lf,%lf]\n",x1,y1);
      for (x2 = x1 + aprox_painter_delta; x2 - b < 1.e-6; x2 += aprox_painter_delta)
        {
          y2 = fabs(f(x2) - bess_aprox.approx_function(x2,n,x,table));
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          //printf("[DEBUG] Point [%lf,%lf]\n",x2,y2);
          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = fabs(f(x2) - bess_aprox.approx_function(x2,n,x,table));
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      painter.setPen ("blue");
//      painter.drawText (0, 40, approx_name);
//      free(x);
//      free(y);
//      free(y_derivatives);
//      free(table);
      break;

  }
  // draw axis
//  painter.setPen (pen_red);
//  painter.drawLine (a, 0, b, 0);
//  painter.drawLine (0, max_y, 0, min_y);

  // restore previously saved Coordinate System
  painter.restore ();

  // render function name
  painter.setPen ("blue");
  painter.drawText (0, 20, f_name);
  //painter.setPen ("blue");
  painter.drawText (0, 40, approx_name);
  painter.drawText (0, 60, QString::number(n));
  painter.drawText (0,80,"max{|F_max|,|F_min|}");
  painter.drawText(160,80,QString::number(fmax(fabs(min_y),fabs(max_y))));
  fflush(stdout);
  delete [] x;
  delete [] y;
  delete [] y_derivatives;
  delete [] table;
}
