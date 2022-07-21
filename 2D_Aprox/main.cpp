#include "mainwindow.h"

#include <QApplication>
#include "openglwidget.h"
#include <QMenuBar>
#include <QAction>
#include <fenv.h>
int main(int argc, char *argv[])
{

    //feenableexcept(FE_UNDERFLOW || FE_OVERFLOW || FE_INVALID || FE_DIVBYZERO);
    int nx = 10;
    int ny = 10;
    double x0 = -1, y0 = -1;
    double x1 = 1, y1 = 1;
    int p = 3;
    int start_func_type;
    char *file_name = new char[512];
    char *buff_line = new char[512];
    double eps = 1e-12;
    FILE* input_file;
    if (argc != 7){
        delete [] file_name;
        delete [] buff_line;
        printf("Not enough arguments. Usage: %s file_name nx ny start_function eps p\n",argv[0]);
        return -1;
    }

    if (   sscanf (argv[1], "%s", file_name) != 1
             || sscanf (argv[2], "%d", &nx) != 1
             || sscanf (argv[3], "%d", &ny) != 1
             || nx <= 0 || ny <= 0
             || sscanf (argv[4], "%d", &start_func_type) != 1
             || (start_func_type > 7 || start_func_type < 0)
             || sscanf (argv[5], "%lf", &eps) != 1
             || sscanf (argv[6], "%d", &p) != 1
           || p < 1) {
        printf("Invalid arguments. Usage: %s file_name nx ny start_function eps p\n",argv[0]);
        return -1;
    }
    if(!(input_file = fopen(file_name,"r"))){
        delete [] file_name;
        delete [] buff_line;
        printf("[ERROR] Unable to open file\n");
        return -1;
    }
    int read_counter = 0;
    while (fgets (buff_line, 512, input_file))
        {
          if (*buff_line == '#')
            continue;

          if (read_counter == 0)
            {
              if (sscanf (buff_line, "%lf %lf", &x0, &y0) == 2)
                 read_counter = 1;
            }
          else
            {
              if (sscanf (buff_line, "%lf %lf", &x1, &y1) == 2)
                {
                 read_counter = 2;
                 break;
                }
            }
        }
    if (read_counter != 2){
        printf("[ERROR] Malformed or corrupted domain config file\n");
        delete [] file_name;
        delete [] buff_line;
        return -1;
    }
    //printf("AAAAAAAAAAAAAAAAAAAAAAAAAA\n");
    fflush(stdout);
    fclose(input_file);
    delete [] file_name;
    delete [] buff_line;
    QApplication a(argc, argv);
    MainWindow* w = new MainWindow();
    QMenuBar* menubar = new QMenuBar();
    QAction* action;
    openglwidget* painter = new openglwidget(w,nx,ny,x0,x1,y0,y1,p,start_func_type,eps);
    action = menubar->addAction("Rotate cw",painter,SLOT( rotate_ccw() ));
    action->setShortcut('8');
    action = menubar->addAction("Rotate ccw",painter,SLOT( rotate_cw() ));
    action->setShortcut('9');
    action = menubar->addAction("Change function",painter,SLOT( change_func() ));
    action->setShortcut('0');
    action = menubar->addAction("Chenge graphs",painter,SLOT( change_graph() ));
    action = menubar->addAction("Double scale",painter,SLOT( double_scale() ));
    action->setShortcut('2');
    action = menubar->addAction("Sub scale",painter,SLOT( sub_scale() ));
    action->setShortcut('3');
//    action = menubar->addAction("Double nx",painter,SLOT( double_nx() ));
//    action = menubar->addAction("Sub nx",painter,SLOT( sub_nx() ));
//    action = menubar->addAction("Double ny",painter,SLOT( double_ny() ));
//    action = menubar->addAction("Sub ny",painter,SLOT( sub_ny() ));
    action = menubar->addAction("Double nx*ny",painter,SLOT( double_nx_ny() ));
    action->setShortcut('4');
    action = menubar->addAction("sub nx*ny",painter,SLOT( sub_nx_ny() ));
    action->setShortcut('5');
    action = menubar->addAction("increment error",painter,SLOT( increment_error() ));
    action->setShortcut('6');
    action = menubar->addAction("decrement error",painter,SLOT( decrement_error() ));
    action->setShortcut('7');
    w->setMenuBar(menubar);
    w->setCentralWidget(painter);
    w->show();
    return a.exec();
}
