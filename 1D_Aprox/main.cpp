
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <fenv.h>
#include "window.h"

int main (int argc, char *argv[])
{
  feenableexcept(FE_UNDERFLOW || FE_OVERFLOW || FE_DIVBYZERO || FE_INVALID);
  QApplication app (argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window);
  QAction *action;

  if (graph_area->parse_command_line (argc, argv))
    {
      QMessageBox::warning (0, "Wrong input arguments!", 
                            "Wrong input arguments!");
      return -1;
    }

  action = tool_bar->addAction ("&Change function", graph_area, SLOT (change_func ()));
  action->setShortcut (QString ("Ctrl+C"));

  action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
  action->setShortcut (QString ("Ctrl+X"));

  action = tool_bar->addAction ("&x2 points", graph_area, SLOT (double_n ()));
  action->setShortcut (QString ("4"));

  action = tool_bar->addAction ("&/2 points", graph_area, SLOT (sub_n ()));
  action->setShortcut (QString ("5"));

  action = tool_bar->addAction ("change aproximation", graph_area, SLOT (change_approx ()));
  action->setShortcut (QString ("1"));

  action = tool_bar->addAction ("Add mistake", graph_area, SLOT (add_mistake ()));
  action->setShortcut (QString ("6"));

  action = tool_bar->addAction ("Sub mistake", graph_area, SLOT (sub_mistake ()));
  action->setShortcut (QString ("7"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");

  window->show ();
  app.exec ();
  delete window;
  return 0;
}
