
#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  int approx_id;
  const char *f_name;
  const char *approx_name;
  double a;
  double b;
  int n;
  double (*f) (double);
  int mistake_scale;
public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);

public slots:
  void change_func ();
  void change_approx();
  void double_n();
  void sub_n();
  void add_mistake();
  void sub_mistake();
protected:
  void paintEvent (QPaintEvent *event);
};

#endif
