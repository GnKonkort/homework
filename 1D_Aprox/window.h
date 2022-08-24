
#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>
#include "besselaproximator.h"
#include "newtonaproximator.h"

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  int mode_id;
  const char *f_name;
  const char *mode_name;
  double a;
  double b;
  int n;
  double (*f) (double);
  BesselAproximator bessel;
  NewtonAproximator newton;
  int mistake_amount;
public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);

public slots:
  void change_func ();
  void change_mode();
  void double_n();
  void sub_n();
  void zoom_in();
  void zoom_out();
  void increase_error();
  void decrease_error();

protected:
  void paintEvent (QPaintEvent *event);
};

#endif
