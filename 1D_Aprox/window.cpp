
#include <QPainter>
#include <stdio.h>

#include "window.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10

static
double f_0 (double x)
{
  (void)x;
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
    return 1/(25*x*x + 1);
}

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

  func_id = 0;
  mode_id = 0;
  mistake_amount = 0;

  change_func ();
  change_mode();
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

  if (   sscanf (argv[1], "%lf", &this->a) != 1
      || sscanf (argv[2], "%lf", &this->b) != 1
      || sscanf (argv[3], "%d", &this->n) != 1
      || sscanf (argv[4], "%d", &func_id) != 1
      || (func_id > 6 || func_id < 0)
      || b - a < 1.e-6
      || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
      || n <= 0)
    return -2;

  return 0;
}

void Window::change_mode(){
    mode_id = (mode_id + 1) % 4;
    switch(mode_id){
        case 0:
            mode_name = "Newton Aproximation";
        break;
        case 1:
            mode_name = "Bessel Aproximation";
        break;
        case 2:
            mode_name = "Real function and both aproximations";
        break;
        case 3:
            mode_name = "Residual for both aproximations";
        break;
    }
    update();
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
      f_name = "f (x) = x*x";
      f = f_2;
      break;
    case 3:
      f_name = "f (x) = x*x*x";
      f = f_3;
      break;
    case 4:
      f_name = "f (x) = x*x*x*x";
      f = f_4;
      break;
    case 5:
      f_name = "f (x) = exp(x)";
      f = f_5;
      break;
    case 6:
      f_name = "f(x) = 1/(25*x*x + 1)";
      f = f_6;
      break;
    }
  update ();
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{  
  double fmin,fmax;
  QRect rec = QApplication::desktop()->screenGeometry();
  QPainter painter (this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  double delta_y, delta_x = (b - a) / n;
  double screen_delta = (b - a)/rec.width();
  QPen pen_black(Qt::black, 0, Qt::SolidLine); 
  QPen pen_red(Qt::red, 0, Qt::SolidLine); 
  QPen pen_blue(Qt::blue, 0, Qt::SolidLine);

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

  // save current Coordinate System
  painter.save ();

  // make Coordinate Transformations
  painter.translate (0.5 * width (), 0.5 * height ());
  painter.scale (width () / (b - a), -height () / (max_y - min_y));
  painter.translate (-0.5 * (a + b), -0.5 * (min_y + max_y));

  // draw approximated line for graph
//  x1 = a;
//  y1 = f (x1);
//  for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
//    {
//      y2 = f (x2);
//      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

//      x1 = x2, y1 = y2;
//    }
//  x2 = b;
//  y2 = f (x2);
//  painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
  switch (mode_id) {
    case(0)://Newton func
        printf("[DEBUG] Newton Aproximation\n");
        painter.setPen(pen_black);
        newton.build_table(n+1,f,a,b,mistake_amount);
        x1 = a;
        y1 = newton.aprox(x1);
        fmax = y1;
        fmin = y1;
        for (x2 = x1 + screen_delta; x2 - b < 1.e-6; x2 += screen_delta)
        {
          y2 = newton.aprox(x2);
          if(fabs(y2) > fmax) fmax = fabs(y2);
          if(fabs(y2) < fmin) fmin = fabs(y2);
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
          x1 = x2, y1 = y2;
        }
        x2 = b;
        y2 = newton.aprox(x2);
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        painter.setPen(pen_black);
          x1 = a;
          y1 = f (x1);
          if(fabs(y1) > fmax) fmax = fabs(y1);
          if(fabs(y1) < fmin) fmin = fabs(y1);
          for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
            {
              y2 = f (x2);
              if(fabs(y2) > fmax) fmax = fabs(y2);
              if(fabs(y2) < fmin) fmin = fabs(y2);
              painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

              x1 = x2, y1 = y2;
            }
          if(fabs(y2) > fmax) fmax = fabs(y2);
          if(fabs(y2) < fmin) fmin = fabs(y2);
          x2 = b;
          y2 = f (x2);
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      break;
    case(1)://Bessel func
      printf("[DEBUG] Bessel Aproximation\n");
      painter.setPen(pen_black);
      bessel.build_table(n+1,a,b,f,mistake_amount);
      x1 = a;
      y1 = bessel.aprox(x1);
      fmax = y1;
      fmin = y1;
      for (x2 = x1 + screen_delta; x2 - b < 1.e-6; x2 += screen_delta)
      {
        y2 = bessel.aprox(x2);
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        x1 = x2, y1 = y2;
      }
      x2 = b;
      y2 = bessel.aprox(x2);
      if(fabs(y2) > fmax) fmax = fabs(y2);
      if(fabs(y2) < fmin) fmin = fabs(y2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      painter.setPen(pen_black);
        x1 = a;
        y1 = f (x1);
        if(fabs(y1) > fmax) fmax = fabs(y1);
        if(fabs(y1) < fmin) fmin = fabs(y1);
        for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
          {
            y2 = f (x2);
            if(fabs(y2) > fmax) fmax = fabs(y2);
            if(fabs(y2) < fmin) fmin = fabs(y2);
            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

            x1 = x2, y1 = y2;
          }
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        x2 = b;
        y2 = f (x2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      break;
    case(2)://Both
      printf("[DEBUG] Both\n");
      painter.setPen(pen_red);
      newton.build_table(n+1,f,a,b,mistake_amount);
      x1 = a;
      y1 = newton.aprox(x1);
      fmax = y1;
      fmin = y1;
      for (x2 = x1 + screen_delta; x2 - b < 1.e-6; x2 += screen_delta)
      {
        y2 = newton.aprox(x2);
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        x1 = x2, y1 = y2;
      }
      x2 = b;
      y2 = newton.aprox(x2);
      if(fabs(y2) > fmax) fmax = fabs(y2);
      if(fabs(y2) < fmin) fmin = fabs(y2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      painter.setPen(pen_blue);
      bessel.build_table(n+1,a,b,f,mistake_amount);
      x1 = a;
      y1 = bessel.aprox(x1);
      if(fabs(y1) > fmax) fmax = fabs(y1);
      if(fabs(y1) < fmin) fmin = fabs(y1);
      for (x2 = x1 + screen_delta; x2 - b < 1.e-6; x2 += screen_delta)
      {
        y2 = bessel.aprox(x2);
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        x1 = x2, y1 = y2;
      }
      x2 = b;
      y2 = bessel.aprox(x2);
      if(fabs(y2) > fmax) fmax = fabs(y2);
      if(fabs(y2) < fmin) fmin = fabs(y2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      painter.setPen(pen_black);
        x1 = a;
        y1 = f (x1);
        if(fabs(y1) > fmax) fmax = fabs(y1);
        if(fabs(y1) < fmin) fmin = fabs(y1);
        for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
          {
            y2 = f (x2);
            if(fabs(y2) > fmax) fmax = fabs(y2);
            if(fabs(y2) < fmin) fmin = fabs(y2);
            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

            x1 = x2, y1 = y2;
          }
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        x2 = b;
        y2 = f (x2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      break;
    case(3)://Residual
      printf("[DEBUG] Residual\n");
      painter.setPen(pen_red);
      newton.build_table(n+1,f,a,b,mistake_amount);
      x1 = a;
      y1 = fabs(f(x1) - newton.aprox(x1));
      fmax = y1;
      fmin = y1;
      for (x2 = x1 + screen_delta; x2 - b < 1.e-6; x2 += screen_delta)
      {
        y2 = fabs(f(x2) - newton.aprox(x2));
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        x1 = x2, y1 = y2;
      }
      x2 = b;
      y2 = fabs(f(x2) - newton.aprox(x2));
      if(fabs(y2) > fmax) fmax = fabs(y2);
      if(fabs(y2) < fmin) fmin = fabs(y2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      painter.setPen(pen_blue);
      bessel.build_table(n+1,a,b,f,mistake_amount);
      x1 = a;
      y1 = fabs(f(x1) - bessel.aprox(x1));
      if(fabs(y1) > fmax) fmax = fabs(y1);
      if(fabs(y1) < fmin) fmin = fabs(y1);
      for (x2 = x1 + screen_delta; x2 - b < 1.e-6; x2 += screen_delta)
      {
        y2 = fabs(f(x2) - bessel.aprox(x2));
        if(fabs(y2) > fmax) fmax = fabs(y2);
        if(fabs(y2) < fmin) fmin = fabs(y2);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        x1 = x2, y1 = y2;
      }
      x2 = b;
      y2 = fabs(f(x2) - bessel.aprox(x2));
      if(fabs(y2) > fmax) fmax = fabs(y2);
      if(fabs(y2) < fmin) fmin = fabs(y2);
      break;

  }
  fflush(stdout);

  // draw axis
  painter.setPen (pen_red);
  painter.drawLine (a, 0, b, 0);
  painter.drawLine (0, max_y, 0, min_y);

  // restore previously saved Coordinate System
  painter.restore ();

  // render function name
  painter.setPen ("blue");
  painter.drawText (0, 20, f_name);
  painter.drawText(0,50,QString("{|F_min|,|F_max|} = {%1,%2}").arg(QString::number(fmin),QString::number(fmax)));
  printf("{|F_min|,|F_max|} = {%.3e,%.3e}\n",fmin,fmax);
  painter.drawText(0,80,QString("n = %1").arg(QString::number(n)));
  printf("n = %d\n",n);
  painter.drawText(0,110,QString("p = %1").arg(QString::number(mistake_amount)));
  painter.drawText(0,140,mode_name);
  painter.drawText(0,170,QString("[a,b] = [%1,%2]").arg(QString::number(a),QString::number(b)));
  printf("p = %d\n",mistake_amount);
  fflush(stdout);
}
void Window::double_n(){
    this->n *= 2;
    update();
}
void Window::sub_n(){
    if(n / 2 <= 2){

    } else {
       this->n /= 2;
    }

    update();
}
void Window::zoom_in(){
    a /= 2;
    b /= 2;
    update();
}
void Window::zoom_out(){
    a *= 2;
    b *= 2;
    update();
}
void Window::increase_error(){
    this->mistake_amount++;
    update();
}
void Window::decrease_error(){
    this->mistake_amount--;
    update();
}
