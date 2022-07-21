#include <stdio.h>
#include <mpich/mpi.h>
#include "libs/MPI_matrix.h"
#include <memory>
#include <fenv.h>
int main(int argc, char** argv){
    feenableexcept(FE_UNDERFLOW || FE_OVERFLOW || FE_INVALID || FE_DIVBYZERO);
    //td::unique_ptr<double[]> A,B,C,backup_matrix,A_row_buff, B_row_buff;
    char filename[128];
    if(argc < 5 || argc > 6){
        printf("Usage: %s n m r s [filename]\n",argv[0]);
        return -1;
    }
    int n,m,r,s,p,k;
    if( !sscanf(argv[1],"%d",&n)||
        !sscanf(argv[2],"%d",&m)||
        !sscanf(argv[3],"%d",&r)||
        !sscanf(argv[4],"%d",&s)){
        printf("Usage: %s n m r s [filename]\n",argv[0]);
        return -1;
    }
    //printf("AAAAAA");
    if(argc == 6){
        if(s == 0){
            if(!sscanf(argv[5],"%s",filename)){
                printf("Usage: %s n m r s [filename]\n",argv[0]);
                return -1;
            }
        }
        else {
            printf("Usage: %s n m r s [filename]\n",argv[0]);
            return -1;
        }
    } 
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&k);
    MPI_Comm_size(MPI_COMM_WORLD,&p);


    int amount_of_rows_in_process = get_max_rows(n,m,p);
    if(k == 0){
        printf("[DEBUG] Total amount of processes: %d\n",p);
        printf("[DEBUG] Total amount of rows: %d\n",n / m);
    }
    printf("[DEBUG][%d] Amount of rows:%d \n",k, amount_of_rows_in_process);

    std::unique_ptr<double[]> A{new double[amount_of_rows_in_process * n * m]},
                            B{new double[amount_of_rows_in_process * n * m]},
                            C{new double[amount_of_rows_in_process * n * m]},
                            backup_matrix{new double[amount_of_rows_in_process * n * m]},
                            A_row_buff{new double[n*m]}, 
                            B_row_buff{new double[n*m]};
    
    if(s != 0){
        if(s == 1){
            init_matrix(A.get(),n,m,p,k,&formula_1);
        } else if(s == 2){
            init_matrix(A.get(),n,m,p,k,&formula_2);
        } else if(s == 3){
            init_matrix(A.get(),n,m,p,k,&formula_3);
        } else if(s == 4){
            init_matrix(A.get(),n,m,p,k,&formula_4);
        }
    } else {
        if(read_matrix(A.get(),n,m,p,k,filename,A_row_buff.get()) == 1){
            printf("[DEBUG] Failed to open file\n");
            return -1;
        }
    }
    memcpy(backup_matrix.get(),A.get(),amount_of_rows_in_process * m * n);
    init_matrix(B.get(),n,m,p,k,&formula_E);
    if(k == 0){
        printf("Initial matrix:\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    print_matrix(A.get(),n,m,p,k,A_row_buff.get(),r);
    MPI_Barrier(MPI_COMM_WORLD);

    memcpy(backup_matrix.get(),A.get(),amount_of_rows_in_process*n*m*8);
    double err = parallel_matrix_norm(A.get(),k,amount_of_rows_in_process*m,n);
    // printf("[DEBUG] Matrix norm = %.3e\n",err);
    double elapsed = MPI_Wtime();
    int status = solve(A.get(),B.get(),n,m,p,k,A_row_buff.get(),B_row_buff.get(),err);
    elapsed = MPI_Wtime() - elapsed;
    if(status == -1){
        printf("[DEBUG] Matrix is irreversable\n");
        return -1;
    }
    // if(k == 0){
    //     printf("Initial matrix after solution:\n");
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // print_matrix(A.get(),n,m,p,k,A_row_buff.get(),r);
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(k == 0){
    //     printf("Reverse matrix after solution:\n");
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // print_matrix(B.get(),n,m,p,k,A_row_buff.get(),r);
    // MPI_Barrier(MPI_COMM_WORLD);

    mpi_matrix_multiplication(B.get(),backup_matrix.get(),C.get(),n,m,p,k);
    if(k == 0){
        printf("Check matrix solution:\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    print_matrix(C.get(),n,m,p,k,A_row_buff.get(),r);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("[DEBUG] Residual = %.3e\n",parallel_matrix_norm(C.get(),k,amount_of_rows_in_process*m,n));
    double residual = parallel_matrix_norm(C.get(),k,amount_of_rows_in_process*m,n);
    if(k == 0){
        printf("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], residual, elapsed, s, n, m);
    }
    MPI_Finalize();
    return 0;
}