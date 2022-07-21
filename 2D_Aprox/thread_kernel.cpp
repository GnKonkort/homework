#include "thread_kernel.h"
#include "aproximator.h"
#include <limits>
void *thread_kernel::thread_function(void *a)
{
    int result = 0;
     //declaring class to gain access to all required operations
    ThreadArgument* args = (ThreadArgument*) a;
    aproximator* aprox = new aproximator(args->common);
    //printf("[THREAD][%d] Succesfully started, entering loop\n",args->k);
//    printf("[DEBUG][%d] args->common->cond = %d\n",args->k,&args->common->cond);
    while(true){
        pthread_mutex_lock(&(args->common->mutex_calculation));
        while(args->common->current_task_status == TaskStatus::NO_TASK || args->common->current_task_status == TaskStatus::TASK_COMPLETE){
            //printf("[DEBUG][%d] Thread is awaiting for new task\n",args->k);
            fflush(stdout);
            pthread_cond_wait(&(args->common->cond),&(args->common->mutex_calculation));
        }
        pthread_mutex_unlock(&(args->common->mutex_calculation));
        if(args->common->current_task_status == TaskStatus::CHANGED_FUNC || args->common->current_task_status == TaskStatus::CHANGED_NXNY){
            int N = (args->common->nx + 1) * (args->common->ny + 1);
            //printf("[DEBUG] (%lf,%lf)<->(%lf,%lf)\n",args->common->x0,args->common->y0,args->common->x1,args->common->y1);
            if(args->k == 0){
                pthread_mutex_lock(&args->common->mutex);
                aprox->allocate_arrays(
                            args->common->nx,
                            args->common->ny,
                            args->p,
                            args->common->I,
                            args->common->A,
                            args->common->B,
                            args->common->U,
                            args->common->X,
                            args->common->R,
                            args->common->V,
                            args->common->buf
                            );
                args->common->LEN = aprox->get_len_msr(args->common->nx,args->common->ny);
                aprox->allocate_msr_matrix(
                            args->common->nx,
                            args->common->ny,
                            args->common->A,
                            args->common->I
                            );
                pthread_mutex_unlock(&args->common->mutex);
            }
            aprox->reduce_sum(args->p,0,1);
            aprox->build_msr_matrix(
                        args->common->nx,
                        args->common->ny,
                        args->common->A,
                        args->common->I,
                        args->p,
                        args->k
                        );
            aprox->allocate_msr_vec(
                        args->common->B,
                        args->common->nx,
                        args->common->ny,
                        args->p,
                        args->k,
                        //args->common->f,
                        args->common->x0,
                        args->common->y0,
                        args->common->x1,
                        args->common->y1
                        );

            aprox->reduce_sum(args->p,0,1);
            if(args->k == 0){
                aprox->normalize_msr(
                            args->common->A,
                            args->common->B,
                            args->common->nx,
                            args->common->ny,
                            args->p,
                            args->k,
                            args->common->x0,
                            args->common->y0,
                            args->common->x1,
                            args->common->y1,
                            args->common->LEN,
                            N
                            );
            }
            aprox->reduce_sum(args->p,0,1);
            result = aprox->solver(
                        (args->common->nx + 1) * (args->common->ny + 1),
                        args->common->A,
                        args->common->I,
                        args->common->B,
                        args->common->X,
                        args->common->R,
                        args->common->U,
                        args->common->V,
                        args->common->eps,
                        1000,
                        args->p,
                        args->k,
                        args->common->buf,
                        &(args->iterations)
                        );
            //printf("[DEBUG] Solving complete!\n");
            fflush(stdout);
            if(result < 0){printf("[THREAD][%d] Unable to solve matrix!\n",args->k);fflush(stdout);}
            aprox->reduce_sum(args->p,0,1);

            pthread_mutex_lock(&args->common->mutex);
            args->thread_status = ThreadStatus::FREE,
            args->common->current_task_status = TaskStatus::TASK_COMPLETE,
            args->common->iterations = args->iterations;
            pthread_mutex_unlock(&args->common->mutex);
            //printf("[DEBUG][%d] Task complete\n",args->k);
            fflush(stdout);
        }
    }
}

void thread_kernel::assign_new_task(int nx, int ny,double x0,double x1, double y0, double y1, double (*f)(double, double), int error_amount,double eps)
{
    pthread_mutex_lock(&common->mutex);
    common->nx = nx;
    common->ny = ny;
    common->f = f;
    common->x0 = x0;
    common->x1 = x1;
    common->y0 = y1;
    common->y1 = y0;
    //common->eps = 1e-14;
    common->eps = eps;
    common->error_amount = error_amount;
    common->f_max = std::numeric_limits<double>::lowest();
    //common->error_amount = error_amount;
    printf("[DEBUG] Assigning new task with error_emount = %d eps = %.3e\n",common->error_amount,common->eps);
    double hx = (x1 - x0) / 64;
    double hy = (y1 - y0) / 64;
    for(int i = 0; i < 64; i++){
        for(int j = 0; j < 64; j++){
            if(f(x0 + i * hx,y0 + j * hy) > common->f_max){
                common->f_max = f(x0 + i * hx,y0 + j * hy);
            }
        }
    }
    //printf("[DEBUG]common->f_max = %lf\n",common->f_max);
    fflush(stdout);
    common->current_task_status = TaskStatus::CHANGED_FUNC;
    printf("[THREAD KERNEL] Assigning new task with parameter nx = %d ny = %d\n",common->nx,common->ny);
    pthread_mutex_unlock(&common->mutex);
    pthread_cond_broadcast(&common->cond);
}

thread_kernel::thread_kernel(int p)
{
    threads = new pthread_t[p];
    arguments = new ThreadArgument[p];
    common = new Databus();
    for(int i = 0; i < p; i++){
        arguments[i].common = common;
        arguments[i].k = i;
        arguments[i].p = p;
        arguments[i].thread_status = ThreadStatus::FREE;
    }
    for(int i = 0; i < p; i++){
        pthread_create(threads + i,0,thread_function, arguments + i);
    }
}
