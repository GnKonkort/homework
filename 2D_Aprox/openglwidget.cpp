#include "openglwidget.h"
#include <QPainter>
#include "math.h"
#include <stdio.h>
#include <unistd.h>
double f0(double x, double y){
    (void)x;
    (void)y;
    return 1;
}
double f1(double x, double y){
    (void)y;
    return x;
}
double f2(double x, double y){
    (void)x;
    return y;
}
double f3(double x, double y){
    return x + y;
}
double f4(double x, double y){
    return sqrt(x*x + y*y);
}
double f5(double x, double y){
    return x*x + y*y;
}
double f6(double x, double y){
    return exp(x*x - y*y);
}
double f7(double x, double y){
    return 1 / (25*(x*x + y*y) + 1);
}


double basis(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y){
    return ((x-x1)*(y3-y2) + (y-y1)*(x2-x3))/((x3-x1)*(y2-y1) - (x2-x1)*(y3-y1)) + 1;
}

int FindX(double x, double a, double b, double nx){
    if (x <= a) return 0;
    if (x >= b) return nx;

    int i = 0;
    int j = nx;
    int z = 0;
    double h = (b - a)/nx;
    while (i != j){
        z = (i+j)/2;
        if (x <= a + z*h)
            j = z;
        else
            i = z+1;
    }

  return i;
}

int FindY(double y, double c, double d, double ny){
    if (y <= c) return 0;
    if (y >= d) return ny;

    int i = 0;
    int j = ny;
    int z = 0;
    double h = (d - c)/ny;
    while (i != j){
        z = (i+j)/2;
        if (y <= c + z*h)
            j = z;
        else
            i = z+1;
    }

    return i-1;
}

double aproximator_function(double x, double y, int nx, int ny, double x0, double y0, double x1, double y1, double* result){
    double dx = (x1 - x0)/nx;
        double dy = (y1 - y0)/ny;
        double tmp1, tmp2;

        int i = FindX(x, x0, x1, nx);
        int j = ny - FindY(y, y0, y1, ny);

        if(fabs(y - y1 + j*dy) < 1.e-14 && fabs(x - x0 - i*dx) < 1.e-14)
            return result[(ny-j)*(nx+1) + i];

        if(fabs(y - y1) < 1.e-14){
            double x1_ = x0 + (i-1)*dx;
            double y1_ = y1 - j*dy;

            double x2_ = x0 + i*dx;
            double y2_ = y1 - j*dy;

            double x3_ = x0 + (i-1)*dx;
            double y3_ = y1 - (j+1)*dy;

            return result[i + (ny-j)*(nx+1)]*basis(x2_, y2_, x1_, y1_, x3_, y3_, x, y) +
                   result[i-1 + (ny-j)*(nx+1)]*basis(x1_, y1_, x2_, y2_, x3_, y3_, x, y);
        }
        else if(fabs(y - y0) < 1.e-14){
            double x1_ = x0 + (i-1)*dx;
            double y1_ = y1 - j*dy;

            double x2_ = x0 + i*dx;
            double y2_ = y1 - j*dy;

            double x3_ = x0 + i*dx;
            double y3_ = y1 - (j-1)*dy;

            return result[i + (ny-j)*(nx+1)]*basis(x2_, y2_, x1_, y1_, x3_, y3_, x, y) +
                   result[i-1 + (ny-j)*(nx+1)]*basis(x1_, y1_, x2_, y2_, x3_, y3_, x, y);
        }
        else if(fabs(x - x0) < 1.e-14){
            double x1_ = x0 + i*dx;
            double y1_ = y1 - (j-1)*dy;

            double x2_ = x0 + i*dx;
            double y2_ = y1 - j*dy;

            double x3_ = x0 + (i+1)*dx;
            double y3_ = y1 - (j-1)*dy;

            return result[i + (ny-j)*(nx+1)]*basis(x2_, y2_, x1_, y1_, x3_, y3_, x, y) +
                   result[i + (ny-(j-1))*(nx+1)]*basis(x1_, y1_, x2_, y2_, x3_, y3_, x, y);
        }
        else if(fabs(x - x1) < 1.e-14){
            double x1_ = x0 + i*dx;
            double y1_ = y1 - (j-1)*dy;

            double x2_ = x0 + i*dx;
            double y2_ = y1 - j*dy;

            double x3_ = x0 + (i-1)*dx;
            double y3_ = y1 - j*dy;

            return result[i + (ny-j)*(nx+1)]*basis(x2_, y2_, x1_, y1_, x3_, y3_, x, y) +
                   result[i + (ny-(j-1))*(nx+1)]*basis(x1_, y1_, x2_, y2_, x3_, y3_, x, y);
        }


        double x1_ = x0 + (i-1)*dx;
        double y1_ = y1 - (j-1)*dy;

        double x2_ = x0 + i*dx;
        double y2_ = y1 - (j-1)*dy;

        double x3_ = x0 + (i-1)*dx;
        double y3_ = y1 - j*dy;

        double x4_ = x0 + i*dx;
        double y4_ = y1 - j*dy;

        tmp1 = pow(x1_-x, 2.) + pow(y1_-y, 2.);
        tmp2 = pow(x4_-x, 2.) + pow(y4_-y, 2.);

        if(tmp1 < tmp2)
            return result[i-1 + (ny-j+1)*(nx+1)]*basis(x1_, y1_, x2_, y2_, x3_, y3_, x, y) +
                   result[i + (ny-j+1)*(nx+1)]*basis(x2_, y2_, x1_, y1_, x3_, y3_, x, y) +
                   result[i-1 + (ny-j)*(nx+1)]*basis(x3_, y3_, x1_, y1_, x2_, y2_, x, y);

        else
            return result[i + (ny-j)*(nx+1)]*basis(x4_, y4_, x2_, y2_, x3_, y3_, x, y) +
                   result[i + (ny-j+1)*(nx+1)]*basis(x2_, y2_, x4_, y4_, x3_, y3_, x, y) +
                   result[i-1 + (ny-j)*(nx+1)]*basis(x3_, y3_, x4_, y4_, x2_, y2_, x, y);
}



//----------------------------------------------------------------------------------------------//


double openglwidget::f_local(double x, double y)
{
    //printf("[DEBUG] f_max = %lf error_amount = %d\n",this->f_max, this->error_amount);
    double a = (this->x1 + this->x0) / (2);
    double b = (this->y1 + this->y0) / (2);
    //printf("[DEBUG] a = %lf b = %lf\n",a, b);
    if(fabs(x - a) < __DBL_EPSILON__ && fabs(y - b) < __DBL_EPSILON__){
        //printf("AAAAAAAAAAAAAAAAAAA");
        printf("[WIDGET] f_max = %lf\n",f_max);
        return this->f_max * 0.1 * this->error_amount + f(x,y);
    }
    return f(x,y);
}

QSize openglwidget::sizeHint()
{
    return QSize(640,640);
}

void openglwidget::initializeGL()
{
    QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
    f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    //glEnable(GL_DEPTH_TEST);
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_CULL_FACE);
    f->glClearColor(0.5,0.5,0.5,1);
}

void openglwidget::paintGL()
{
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    QPainter p(this);
    QFont font("Times New Roman", 15);

    glLoadIdentity();

    glRotatef(angle_psi,1,0,0);
    glRotatef(angle_fi,0,0,1);
    glScalef(scale,scale,scale);
//    QPainter painter(this);
//    //painter.beginNativePainting();
//    //glClear(GL_COLOR_BUFFER_BIT);
//    painter.setPen(Qt::black);
//    painter.setFont(QFont("Arial", 16));
//    painter.drawText(0, 0, width(), height(), Qt::AlignCenter, "Hello World!");
//    //painter.endNativePainting();


    //this->fill_grid();
    glBegin(GL_TRIANGLES);

    for(int i = 0; i < 64 * 64 * 2; i++){



        glVertex3f(painting_grid->triangles[i].a.x,painting_grid->triangles[i].a.y,painting_grid->triangles[i].a.z);
        glVertex3f(painting_grid->triangles[i].b.x,painting_grid->triangles[i].b.y,painting_grid->triangles[i].b.z);
        glVertex3f(painting_grid->triangles[i].c.x,painting_grid->triangles[i].c.y,painting_grid->triangles[i].c.z);
    }

    glEnd();
    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3f(0,0,0);
    glVertex3f(10,0,0);
    glColor3f(0,1,0);
    glVertex3f(0,0,0);
    glVertex3f(0,10,0);
    glColor3f(0,0,1);
    glVertex3f(0,0,0);
    glVertex3f(0,0,10);
    glColor3f(1,1,1);
    glEnd();
    p.setPen(Qt::red);
    p.setFont(font);
    p.drawText(30,30,QString("nx:"));
    p.drawText(60,30,QString::number(this->nx));
    p.drawText(30,50,QString("ny:"));
    p.drawText(60,50,QString::number(this->ny));
    p.drawText(30,70,QString("p:"));
    p.drawText(60,70,QString::number(this->error_amount));
    p.drawText(30,90,QString("{|F_max|,|F_min|}:"));
    p.drawText(200,90,QString("{%1,%2}").arg(QString::number(this->f_max),QString::number(this->f_min)));
    p.drawText(30,110,QString("Function ID: %1").arg(QString::number(this->func_mode)));
//    QPainter painter(this);
//    painter.beginNativePainting();
//    glClear(GL_COLOR_BUFFER_BIT);
//    painter.endNativePainting();

//    painter.setPen(Qt::black);
//    painter.setFont(QFont("Arial", 16));
//    painter.drawText(0, 0, width(), height(), Qt::AlignCenter, "Hello World!");
}

openglwidget::openglwidget(QWidget* parent, int nx, int ny, double x0, double x1, double y0, double y1, int p, int start_function,double eps)
{
    printf("[PAINTER] Launching painter...\n");
    fflush(stdout);
    this->thread_overwatch = new QTimer(this);
    this->eps = eps;
    connect(thread_overwatch, &QTimer::timeout, this, &openglwidget::check_threads);
    thread_overwatch->start(1000);
    this->painting_grid = new Grid(64,64,x0,y0,x1,y1);
    //this->f = f0;
    func_mode = (FuncMode)start_function;
    switch (func_mode) {
    case 0:
        this->f = f0;
        break;
    case 1:
        this->f = f1;
        break;
    case 2:
        this->f = f2;
        break;
    case 3:
        this->f = f3;
        break;
    case 4:
        this->f = f4;
        break;
    case 5:
        this->f = f5;
        break;
    case 6:
        this->f = f6;
        break;
    case 7:
        this->f = f7;
        break;
    default:
        break;
    }
    this->graph_mode = GraphMode::REAL_FUNC;
    //this->func_mode = FuncMode::CONST;
    this->angle_fi = 45;
    this->angle_psi = -45;
    this->scale = 0.5;
    this->nx = nx;
    this->ny = ny;
    this->x0 = x0;
    this->x1 = x1;
    this->y0 = y0;
    this->y1 = y1;
    this->p = p;
    this->error_amount = 0;
    this->thread_controller = new thread_kernel(p);
    this->X_local = new double[(nx + 1)*(ny + 1)];
    this->changes_blocked = true;
    this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    this->fill_grid();
}

void openglwidget::fill_grid()
{
    this->f_min = std::numeric_limits<double>::max();
    this->f_max = std::numeric_limits<double>::min();
    double hx = (x1 - x0) / 64;
    double hy = (y1 - y0) / 64;
    for(int i = 0; i < 64; i++){
        for(int j = 0; j < 64; j++){
            if(this->f(x0 + i * hx,y0 + j * hy) > this->f_max){
                this->f_max = f(x0 + i * hx,y0 + j * hy);
            }
        }
    }
    switch (this->graph_mode) {
    case GraphMode::REAL_FUNC:
        printf("[DEBUG] Painting real func\n");
        fflush(stdout);
        for(int i = 0; i < 2 * 64 * 64; i++){
            this->painting_grid->triangles[i].a.z = f_local(this->painting_grid->triangles[i].a.x,this->painting_grid->triangles[i].a.y);
            this->painting_grid->triangles[i].b.z = f_local(this->painting_grid->triangles[i].b.x,this->painting_grid->triangles[i].b.y);
            this->painting_grid->triangles[i].c.z = f_local(this->painting_grid->triangles[i].c.x,this->painting_grid->triangles[i].c.y);
        }
        this->residual = 0;
        break;
    case GraphMode::APROX_FUNC:
        printf("[DEBUG] Painting aprox func\n");\
        //printf("[DEBUG] Aproximation id:%d\n",this->func_mode);
        fflush(stdout);
        for(int i = 0; i < 2 * 64 * 64; i++){
            this->painting_grid->triangles[i].a.z = aproximator_function(
                        this->painting_grid->triangles[i].a.x,
                        this->painting_grid->triangles[i].a.y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        );
            this->painting_grid->triangles[i].b.z = aproximator_function(
                        this->painting_grid->triangles[i].b.x,
                        this->painting_grid->triangles[i].b.y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        );
            this->painting_grid->triangles[i].c.z = aproximator_function(
                        this->painting_grid->triangles[i].c.x,
                        this->painting_grid->triangles[i].c.y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        );
        }
        this->residual = 0;
        break;
    case GraphMode::RESIDUAL:
        //double residual;
        this->residual = 0;
        for(int i = 0; i < 2 * 64 * 64; i++){
            if(fabs(aproximator_function(
                        this->painting_grid->triangles[i].center_of_mass().x,
                        this->painting_grid->triangles[i].center_of_mass().y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        ) - f_local(this->painting_grid->triangles[i].center_of_mass().x,
                              this->painting_grid->triangles[i].center_of_mass().y)) > residual){
                    this->residual = fabs(aproximator_function(
                                        this->painting_grid->triangles[i].center_of_mass().x,
                                        this->painting_grid->triangles[i].center_of_mass().y,
                                        this->nx,
                                        this->ny,
                                        this->x0,
                                        this->y0,
                                        this->x1,
                                        this->y1,
                                        this->X_local
                                        ) - f_local(this->painting_grid->triangles[i].center_of_mass().x,
                                              this->painting_grid->triangles[i].center_of_mass().y));
            }
        }
        printf("[DEBUG] Painting residual func\n");
        printf("Residual: %.3e\n",residual);
        fflush(stdout);
        for(int i = 0; i < 2 * 64 * 64; i++){
            this->painting_grid->triangles[i].a.z = fabs(aproximator_function(
                        this->painting_grid->triangles[i].a.x,
                        this->painting_grid->triangles[i].a.y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        ) - f_local(this->painting_grid->triangles[i].a.x,
                              this->painting_grid->triangles[i].a.y));
            this->painting_grid->triangles[i].b.z = fabs(aproximator_function(
                        this->painting_grid->triangles[i].b.x,
                        this->painting_grid->triangles[i].b.y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        ) - f_local(this->painting_grid->triangles[i].b.x,
                              this->painting_grid->triangles[i].b.y));
            this->painting_grid->triangles[i].c.z = fabs(aproximator_function(
                        this->painting_grid->triangles[i].c.x,
                        this->painting_grid->triangles[i].c.y,
                        this->nx,
                        this->ny,
                        this->x0,
                        this->y0,
                        this->x1,
                        this->y1,
                        this->X_local
                        ) - f_local(this->painting_grid->triangles[i].c.x,
                              this->painting_grid->triangles[i].c.y));
        }
        break;
    default:
        break;
    }
    this->f_min = std::numeric_limits<double>::max();
    this->f_max = std::numeric_limits<double>::min();
    for(int i = 0; i < 2 * 64 * 64; i++){
        if(fabs(painting_grid->triangles[i].a.z) > f_max){
            f_max = fabs(painting_grid->triangles[i].a.z);
        }
        if(fabs(painting_grid->triangles[i].b.z) > f_max){
            f_max = fabs(painting_grid->triangles[i].b.z);
        }
        if(fabs(painting_grid->triangles[i].c.z) > f_max){
            f_max = fabs(painting_grid->triangles[i].c.z);
        }

        if(fabs(painting_grid->triangles[i].a.z) < f_min){
            f_min = fabs(painting_grid->triangles[i].a.z);
        }
        if(fabs(painting_grid->triangles[i].b.z) < f_min){
            f_min = fabs(painting_grid->triangles[i].b.z);
        }
        if(fabs(painting_grid->triangles[i].c.z) < f_min){
            f_min = fabs(painting_grid->triangles[i].c.z);
        }
    }
    printf("{|F_min|,|F_max| = {%.3e,%.3e}\n",f_min,f_max);
}

void openglwidget::rotate_cw()
{
    this->angle_fi += 15;
    update();
}

void openglwidget::rotate_ccw()
{
    this->angle_fi -=15;
    update();
}

void openglwidget::double_scale()
{
    this->scale *= 2;
    update();
}

void openglwidget::sub_scale()
{
    this->scale /= 2;
    update();
}

void openglwidget::double_nx()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
        this->changes_blocked = true;
        this->nx *= 2;
        //printf("[DEBUG] nx = %d\n",nx);
        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
//    update();
}

void openglwidget::sub_nx()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    if(nx > 2){
        this->changes_blocked = true;
        this->nx /= 2;
        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    }
//    update();
}

void openglwidget::double_ny()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
        this->changes_blocked = true;
        this->ny *= 2;
        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    //    update();
}

void openglwidget::double_nx_ny()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    this->changes_blocked = true;
    this->ny *= 2;
    this->nx *= 2;
    this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
}
void openglwidget::sub_nx_ny()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    this->changes_blocked = true;
    if(ny > 2){
        this->ny /= 2;
    }
    if(nx > 2){
        this->nx /= 2;
    }
    this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
}

void openglwidget::sub_ny()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    if(ny >= 2){
        this->changes_blocked = true;
        this->ny /= 2;
        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    }
//    update();
}

void openglwidget::change_func()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    this->func_mode = (FuncMode)((this->func_mode + 1) % 8);
    switch (func_mode) {
    case 0:
        this->f = f0;
        break;
    case 1:
        this->f = f1;
        break;
    case 2:
        this->f = f2;
        break;
    case 3:
        this->f = f3;
        break;
    case 4:
        this->f = f4;
        break;
    case 5:
        this->f = f5;
        break;
    case 6:
        this->f = f6;
        break;
    case 7:
        this->f = f7;
        break;
    default:
        break;
    }

//    if(graph_mode == GraphMode::APROX_FUNC || graph_mode == GraphMode::RESIDUAL){
//         this->changes_blocked = true;
//        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
//    } else {
//        this->fill_grid();
//    }
    this->changes_blocked = true;
    this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    update();
}

void openglwidget::change_graph()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    this->graph_mode = (GraphMode)((this->graph_mode + 1) % 3);
//    if(graph_mode == GraphMode::APROX_FUNC || graph_mode == GraphMode::RESIDUAL){
//        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
//    } else {
//        this->fill_grid();
//    }
    this->fill_grid();
    update();
}

void openglwidget::check_threads()
{
    if(thread_controller->common->current_task_status == TaskStatus::TASK_COMPLETE){
        printf("[THREAD OVERWATCH] Task completed succesfully, retrieving results\n");
        thread_controller->common->current_task_status = TaskStatus::NO_TASK;
        delete [] this->X_local;
        X_local = new double [(nx + 1)*(ny + 1)];
        memcpy(X_local,thread_controller->common->X,8 * (nx + 1)*(ny + 1));
        fill_grid();
        this->changes_blocked = false;
    }
    update();
}

void openglwidget::increment_error()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }
    this->error_amount++;
    /*if(graph_mode == GraphMode::APROX_FUNC || graph_mode == GraphMode::RESIDUAL){
        this->changes_blocked = true;
        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    } else {
        this->fill_grid();
    }*/
    this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    //update();
}

void openglwidget::decrement_error()
{
    if(this->changes_blocked == true){
        printf("[THREAD KERNEL] Threads are currently busy, no changes will be applied\n");
        return;
    }

    this->error_amount--;
    if(graph_mode == GraphMode::APROX_FUNC || graph_mode == GraphMode::RESIDUAL){
        this->changes_blocked = true;
        this->thread_controller->assign_new_task(nx,ny,x0,x1,y0,y1,f,error_amount,this->eps);
    } else {
        this->fill_grid();
    }
    update();
}
