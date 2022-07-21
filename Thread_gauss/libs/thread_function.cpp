#include "thread_function.h"
#include "matrix.h"
#include <limits>
// int solve(double *A, double *B, int n, int m, int k, int p, pthread_mutex_t* mutex,pthread_barrier_t* barrier){

// }

double get_time() {
  struct rusage r;
  getrusage(RUSAGE_THREAD, &r);

  return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1.e6;
}

double get_full_time() {
  struct timeval r;
  gettimeofday(&r, 0);

  return r.tv_sec + r.tv_usec / 1000000.0;
}


void *thread_function(void *g){
    args* a = (args*) g;
    cpu_set_t cpu;
    int nprocs = get_nprocs();
    CPU_ZERO(&cpu);
    CPU_SET(nprocs - 1 - (a->k), &cpu);
    sched_setaffinity(getpid(), sizeof(cpu), &cpu);
    printf("[DEBUG] Thread %d started\n",a->k);

    

    int amount_of_blocks = a->n / a->m;
    int l = a->n % a->m;
    int n = a->n;
    int m = a->m;
    int k = a->k;
    int p = a->p;
    std::unique_ptr<double[]> block1{new double[m*m]};
    std::unique_ptr<double[]> block2{new double[m*m]};
    std::unique_ptr<double[]> block3{new double[m*m]};

    std::unique_ptr<double[]> row_buffer{new double[n*m]};
    std::unique_ptr<double[]> row_buffer2{new double[n*m]};
    a->elapsed_thread = get_time();
    a->wall_clock = get_full_time();
    for(int i = 0; i < amount_of_blocks; i++){
        double buffer;
        for(int j = k; j < amount_of_blocks; j+=p){
            if(j < i) continue;
            get_block(a->A,n,m,j,i,block1.get());
            //printf("[DEBUG] getting block (%d,%d)\n",j,i);
            if(inverse_matrix(block1.get(),block2.get(),m,a->norm) == 0){
                if((buffer = norma(block2.get(),m,m)) < a->reverse_block_norm){
                    a->reverse_block_norm = buffer;
                    a->reverse_block_pos = j;
                    //put_block(a->reverse_block_buffer,m,m,0,0,block2.get());
                }
            }
        }
        pthread_barrier_wait(a->barrier);
        //Leader-thread should now find, which block to pick
        if(k == i % p){
            for(int j = 0; j < p; j++){
                if(a->a[j].reverse_block_norm < a->reverse_block_norm){
                    a->reverse_block_norm = a->a[j].reverse_block_norm;
                    a->reverse_block_pos = a->a[j].reverse_block_pos;
                }
            }
            for(int j = 0; j < p; j++){
                a->a[j].reverse_block_norm = a->reverse_block_norm;
                a->a[j].reverse_block_pos = a->reverse_block_pos;
            }
        }
        pthread_barrier_wait(a->barrier);
        printf("[DEBUG] reverse_block_pos = %d reverse_block_norm = %.3e\n",a->reverse_block_pos,a->reverse_block_norm);
        if(a->reverse_block_pos == -1){
            a->error = -1;
            return 0;
        }
        //Now all threads are aware, which row is main. Now they should switch rows
        for(int j = k; j < amount_of_blocks; j+=p){
                get_block(a->A,n,m,i,(j) % amount_of_blocks,block1.get());
                get_block(a->A,n,m,a->reverse_block_pos,(j) % amount_of_blocks,block2.get());
                put_block(a->A,n,m,i,(j) % amount_of_blocks,block2.get());
                put_block(a->A,n,m,a->reverse_block_pos,(j) % amount_of_blocks,block1.get());

                get_block(a->B,n,m,i,(j) % amount_of_blocks,block1.get());
                get_block(a->B,n,m,a->reverse_block_pos,(j) % amount_of_blocks,block2.get());
                put_block(a->B,n,m,i,(j) % amount_of_blocks,block2.get());
                put_block(a->B,n,m,a->reverse_block_pos,(j) % amount_of_blocks,block1.get());
        }
        if(l != 0 && k == amount_of_blocks % p){
            get_block(a->A,n,m,i,amount_of_blocks,block1.get());
            get_block(a->A,n,m,a->reverse_block_pos,amount_of_blocks,block2.get());
            put_block(a->A,n,m,i,amount_of_blocks,block2.get());
            put_block(a->A,n,m,a->reverse_block_pos,amount_of_blocks,block1.get());

            get_block(a->B,n,m,i,amount_of_blocks,block1.get());
            get_block(a->B,n,m,a->reverse_block_pos,amount_of_blocks,block2.get());
            put_block(a->B,n,m,i,amount_of_blocks,block2.get());
            put_block(a->B,n,m,a->reverse_block_pos,amount_of_blocks,block1.get());
        }
        pthread_barrier_wait(a->barrier);
        // if(k == 0){
        //     printf("[DEBUG] Matrix A after switching rows %d and %d\n",i,a->reverse_block_pos);
        //     print_matrix(a->A, n, n, n);
        //     printf("[DEBUG] Matrix B after switching rows %d and %d\n",i,a->reverse_block_pos);
        //     print_matrix(a->B, n, n, n);
        // }
        if(k == i % p){
            get_block(a->A,n,m,i,i,block1.get());
            if(inverse_matrix(block1.get(),block2.get(),m,a->norm) == 0){
                put_block(a->reverse_block_buffer,m,m,0,0,block2.get());
                a->error = 0;
            } else {
                a->error = -1;
            }
        }
        pthread_barrier_wait(a->barrier);
        //pthread_mutex_lock(a->mutex);
        if(a->a[i % p].error == -1){
            a->error = -1;
            return 0;
        }
        get_block(a->reverse_block_buffer,m,m,0,0,block2.get());
        //pthread_mutex_unlock(a->mutex);
        // printf("[DEBUG] Reverse block:\n");
        // print_matrix(block2.get(),m,m,10);

        //Reverse block recieved, now it is time to normalize main row
        if(i % p == k){
            ones(block1.get(),m);
            put_block(a->A,n,m,i,i,block1.get());
        }
        for(int j = i + k + 1; j < amount_of_blocks; j+=p){
            get_block(a->A, n, m, i, j, block1.get());
            multiplication(block2.get(),block1.get(),block3.get(), m, m, m);
            put_block(a->A, n, m, i, j, block3.get());
        }
        if(l != 0 && k == amount_of_blocks % p){
            //printf("[DEBUG][%d] Normalizing border block of row\n",k);
            get_block(a->A, n, m, i, amount_of_blocks, block1.get());
            multiplication(block2.get(),block1.get(),block3.get(), m, m, l);
            put_block(a->A, n, m, i, amount_of_blocks, block3.get());
        }
        for(int j = 0 + k; j < i; j += p){
            get_block(a->B, n, m, i, j, block1.get());
            multiplication(block2.get(),block1.get(),block3.get(), m, m, m);
            put_block(a->B, n, m, i, j, block3.get());
        }
        if(a->reverse_block_pos % p == k){
            get_block(a->B, n, m, i, a->reverse_block_pos, block1.get());
            multiplication(block2.get(),block1.get(),block3.get(), m, m, m);
            put_block(a->B, n, m, i, a->reverse_block_pos, block3.get());
        }
        pthread_barrier_wait(a->barrier);

        //Row is normalized, now it is time for every thread to copy normalized row to their buffer
        for(int j = 0; j < amount_of_blocks; j++){
            get_block(a->A,n,m,i,(j + k) % amount_of_blocks,block1.get());
            put_block(row_buffer.get(),n,m,0,(j + k) % amount_of_blocks,block1.get());

            get_block(a->B,n,m,i,(j + k) % amount_of_blocks,block1.get());
            put_block(row_buffer2.get(),n,m,0,(j + k) % amount_of_blocks,block1.get());
            if((j + k) % amount_of_blocks == 0 && l != 0){
                get_block(a->A,n,m,i,amount_of_blocks,block1.get());
                put_block(row_buffer.get(),n,m,0,amount_of_blocks,block1.get());

                get_block(a->B,n,m,i,amount_of_blocks,block1.get());
                put_block(row_buffer2.get(),n,m,0,amount_of_blocks,block1.get());
            }
        }
        //All threads now have the same row, the can independently start subtracting one row from another
        // pthread_mutex_lock(a->mutex);
        // printf("[DEBUG] Got A matrix buffer-row:\n");
        // print_matrix(row_buffer.get(),m,n,n);
        // printf("[DEBUG] Got B matrix buffer-row:\n");
        // print_matrix(row_buffer2.get(),m,n,n);
        // pthread_mutex_unlock(a->mutex);
        for(int j = i + 1; j < amount_of_blocks; j++){
            if(k == j % p){
                get_block(a->A,n,m,j,i,block1.get());
                zero(block2.get(),m,m);
                put_block(a->A,n,m,j,i,block2.get());
                for(int z = i+1; z < amount_of_blocks; z++){
                    get_block(row_buffer.get(),n,m,0,z,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                    get_block(a->A,n,m,j,z,block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(a->A,n,m,j,z,block2.get());
                }
                if(l != 0){
                    get_block(row_buffer.get(),n,m,0,amount_of_blocks,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,l);
                    get_block(a->A,n,m,j,amount_of_blocks,block2.get());
                    substraction(block2.get(),block3.get(),m,l);
                    put_block(a->A,n,m,j,amount_of_blocks,block2.get());
                }
                for(int z = 0; z < i; z++){
                    get_block(row_buffer2.get(),n,m,0,z,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                    get_block(a->B,n,m,j,z,block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(a->B,n,m,j,z,block2.get());
                }
                get_block(row_buffer2.get(),n,m,0,a->reverse_block_pos,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                get_block(a->B,n,m,j,a->reverse_block_pos,block2.get());
                substraction(block2.get(),block3.get(),m,m);
                put_block(a->B,n,m,j,a->reverse_block_pos,block2.get());
            }
        }
        if(l != 0 && k == amount_of_blocks % p){
            get_block(a->A,n,m,amount_of_blocks,i,block1.get());
            zero(block2.get(),l,m);
            put_block(a->A,n,m,amount_of_blocks,i,block2.get());

            for(int z = i+1; z < amount_of_blocks; z++){
                get_block(row_buffer.get(),n,m,0,z,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),l,m,m);
                get_block(a->A,n,m,amount_of_blocks,z,block2.get());
                substraction(block2.get(),block3.get(),l,m);
                put_block(a->A,n,m,amount_of_blocks,z,block2.get());
            }
            if(l != 0){
                get_block(row_buffer.get(),n,m,0,amount_of_blocks,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),l,m,l);
                get_block(a->A,n,m,amount_of_blocks,amount_of_blocks,block2.get());
                substraction(block2.get(),block3.get(),l,l);
                put_block(a->A,n,m,amount_of_blocks,amount_of_blocks,block2.get());
            }
            for(int z = 0; z < i; z++){
                get_block(row_buffer2.get(),n,m,0,z,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),l,m,m);
                get_block(a->B,n,m,amount_of_blocks,z,block2.get());
                substraction(block2.get(),block3.get(),l,m);
                put_block(a->B,n,m,amount_of_blocks,z,block2.get());
            }
            get_block(row_buffer2.get(),n,m,0,a->reverse_block_pos,block2.get());
            multiplication(block1.get(),block2.get(),block3.get(),l,m,m);
            get_block(a->B,n,m,amount_of_blocks,a->reverse_block_pos,block2.get());
            substraction(block2.get(),block3.get(),m,m);
            put_block(a->B,n,m,amount_of_blocks,a->reverse_block_pos,block2.get());
            
        }
        a->reverse_block_norm = std::numeric_limits<double>::max();
        a->reverse_block_pos = -1;
        pthread_barrier_wait(a->barrier);
        // if(k == 0){
        //     printf("[DEBUG] Matrix A after step %d\n",i);
        //     print_matrix(a->A,n,n,10);
        //     printf("[DEBUG] Matrix B after step %d\n",i);
        //     print_matrix(a->B,n,n,10);
        // }
    }
    if(l != 0){
        if(amount_of_blocks % p == k){
            get_block(a->A,n,m,amount_of_blocks,amount_of_blocks,block1.get());
            if(inverse_matrix(block1.get(),block2.get(),l,a->norm) == 0){
                //put_block_raw(a->reverse_block_buffer,n,m,0,0,block2.get(),l,l);
                put_block(a->reverse_block_buffer,l,l,0,0,block2.get());
                a->error = 0;
            } else {
                a->error = -1;
            }
        }
        pthread_barrier_wait(a->barrier);
        if(a->a[amount_of_blocks % p].error == -1){a->error = -1; return 0;}
        get_block_raw(a->reverse_block_buffer,l,l,0,0,block2.get(),l,l);
        // pthread_mutex_lock(a->mutex);
        // printf("[DEBUG][%d] Got corner reverse block:\n",k);
        // print_matrix(block2.get(),l,l,l);
        // pthread_mutex_unlock(a->mutex);
        if(amount_of_blocks % p == k){
             ones(block1.get(),l);
             put_block_raw(a->A,n,m,amount_of_blocks,amount_of_blocks,block1.get(),l,l);
        }
        for(int j = 0 + k; j < amount_of_blocks; j+=p){
            get_block(a->B, n, m, amount_of_blocks, j, block1.get());
            multiplication(block2.get(),block1.get(),block3.get(), l, l, m);
            put_block(a->B, n, m, amount_of_blocks, j, block3.get());
        }
        if(k == amount_of_blocks % p){
            get_block(a->B, n, m, amount_of_blocks, amount_of_blocks, block1.get());
            multiplication(block2.get(),block1.get(),block3.get(), l, l, l);
            put_block(a->B, n, m, amount_of_blocks, amount_of_blocks, block3.get());
        }
    }
    pthread_barrier_wait(a->barrier);
    //Reverse step, ignoring matrix A, writing only to matrix B.
    // if(k == 0){
    //     printf("[DEBUG] Matrix A after forward stage:\n");
    //     print_matrix(a->A,n,n,10);
    //     printf("[DEBUG] Matrix B after forward stage:\n");
    //     print_matrix(a->B,n,n,10);
    // }
    if(l != 0){
        for(int j = 0; j < amount_of_blocks; j++){
            get_block_raw(a->B,n,m,amount_of_blocks,(j + k) % amount_of_blocks, block1.get(),l,m);
            put_block_raw(row_buffer2.get(),n,m,0,(j+k) % amount_of_blocks,block1.get(),l,m);
            if((j + k) % amount_of_blocks == 0){
                get_block_raw(a->B,n,m,amount_of_blocks,amount_of_blocks, block1.get(),l,l);
                put_block_raw(row_buffer2.get(),n,m,0,amount_of_blocks,block1.get(),l,l);
            }
        }
        // printf("[DEBUG] Got bottom row buffer\n");
        // print_matrix(row_buffer2.get(),l,n,10);

        for(int j = 0; j < amount_of_blocks; j++){
            if(k == j % p){
                get_block(a->A,n,m,j,amount_of_blocks,block1.get());
                for(int z = 0; z < amount_of_blocks; z++){
                    get_block_raw(row_buffer2.get(),n,m,0,z,block2.get(),l,m);
                    multiplication(block1.get(),block2.get(),block3.get(),m,l,m);
                    get_block(a->B,n,m,j,z,block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(a->B,n,m,j,z,block2.get());
                }
                get_block_raw(row_buffer2.get(),n,m,0,amount_of_blocks,block2.get(),m,l);
                multiplication(block1.get(),block2.get(),block3.get(),m,l,l);
                get_block(a->B,n,m,j,amount_of_blocks,block2.get());
                substraction(block2.get(),block3.get(),m,l);
                put_block(a->B,n,m,j,amount_of_blocks,block2.get());
            }
        }

    }
    pthread_barrier_wait(a->barrier);
    // if(k == 0){
    //     printf("[DEBUG] Matrix A after bottom reverse stage:\n");
    //     print_matrix(a->A,n,n,10);
    //     printf("[DEBUG] Matrix B after bottom reverse stage:\n");
    //     print_matrix(a->B,n,n,10);
    // }
    for(int i = amount_of_blocks - 1; i >= 0; i--){
        for(int j = 0; j < amount_of_blocks; j++){
            get_block(a->B,n,m,i,(j + k) % amount_of_blocks, block1.get());
            put_block(row_buffer2.get(),n,m,0,(j+k) % amount_of_blocks,block1.get());
            if((j + k) % amount_of_blocks == 0){
                get_block(a->B,n,m,i,amount_of_blocks, block1.get());
                put_block(row_buffer2.get(),n,m,0,amount_of_blocks,block1.get());
            }
        }
        // printf("[DEBUG] Got row buffer\n");
        // print_matrix(row_buffer2.get(),m,n,10);
        for(int j = 0; j < i; j++){
            if(k == j % p){
                get_block(a->A,n,m,j,i,block1.get());
                for(int z = 0; z < amount_of_blocks; z++){
                    get_block(row_buffer2.get(),n,m,0,z,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                    get_block(a->B,n,m,j,z,block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(a->B,n,m,j,z,block2.get());
                }
                if(l != 0){
                    get_block(row_buffer2.get(),n,m,0,amount_of_blocks,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,l);
                    get_block(a->B,n,m,j,amount_of_blocks,block2.get());
                    substraction(block2.get(),block3.get(),m,l);
                    put_block(a->B,n,m,j,amount_of_blocks,block2.get());
                }
            }
        }
        pthread_barrier_wait(a->barrier);
    }
    a->wall_clock = get_full_time() - a->wall_clock;
    a->elapsed_thread = get_time() - a->elapsed_thread;
    return 0;
}