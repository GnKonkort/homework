#include "MPI_matrix.h"
#include "matrix.h"
int max(int a, int b)
{
    return (a > b ? a : b);
}

int min(int a, int b)
{
    return (a < b ? a : b);
}

double formula_1(int i, int j, int n){
    return n - max(i, j);
}
double formula_2(int i, int j, int n){
    return max(i + 1, j + 1);
}
double formula_3(int i, int j, int n){
    return fabs(i-j);
}
double formula_4(int i, int j, int n){
    return 1.0/(i+j+1);
}
double formula_E(int i, int j, int n){
    return (double)(i==j);
}

int get_max_rows(int n, int m, int p){
    int blocks = n / m + (n % m == 0 ? 0 : 1);
    return blocks / p + (blocks % p == 0 ? 0 : 1);
}
int row_global2local(int n, int m, int p, int k, int i_glob){
    int block_row_glob = i_glob / m;
    int block_row_loc = block_row_glob / p;
    return block_row_loc*m + i_glob % m;
}
int row_local2global(int n, int m, int p, int k, int i_loc){
    int block_row_loc = i_loc / m;
    return k * m + p * m * block_row_loc + i_loc % m;
}
int get_rows_block(int n, int m, int p, int k)
{
    int blocks = (n + m - 1) / m;

    return ((blocks % p == 0 || blocks % p < k) ? blocks / p : blocks / p + 1); // (!!исправил!!)
}
int get_rows(int n, int m, int p, int k){
    int blocks = (n + m - 1)/m;
    return ((blocks % p == 0)||(blocks % p > k) ? blocks / p : blocks / p + 1);
}
double parallel_matrix_norm(double* a,int k, int rows, int n){
    double norm = 0;
    double buf = 0;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    //double* norms;
    //norms = (double*)calloc(comm_size,sizeof(double));
    std::unique_ptr<double[]> norms{new double[comm_size]}; 
    for(int i = 0; i < rows; i++){
        buf = 0;
        for(int j = 0; j < n; j++){
            buf += fabs(a[i*n + j]);
        }
        if(buf > norm){
            norm = buf;
        }
    }
    MPI_Gather(&norm,1,MPI_DOUBLE,norms.get(),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    if(k == 0){
        for(int i = 0; i < comm_size; i++){
            if(norms[i] > norm){
                norm = norms[i];
            }
        }
    }
    MPI_Bcast(&norm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    return norm;
}
int read_array(FILE* fp,double* a, int len){
    int i;
    for(i = 0; i < len; i++){
        if(fscanf(fp,"%lf",a+i) != 1){
            return -2;
        }
    }
    return 0;
}

int print_array(const double *a, int n, int m, int printed, int max_print)
{
    int i, j;
    int can_print = max_print - printed;
    int rows = min(m, can_print);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; ((j < max_print)&&(j < n)); j++)
        {
            printf("%lf ", a[i * n + j]);
        }

        printf("\n");
    }

    return rows;
}

void print_matrix(const double *a, int n, int m, int p, int k, double *buf, int max_print)
{
    int main_k = 0;
    int block;
    int max_blocks = (n + m - 1) / m;
    int printed = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    int block_loc = 0;
    int owner;
    //int blocks_in_process = get_rows_block(n, m, p, k);

    for (block = 0; block < max_blocks; block++)
    {
        owner = block % p;
        block_loc = block / p;
        // if(k == 0){
        //     printf("[DEBUG] block_loc = %d\n[DEBUG] block_owner = %d\n",block_loc,owner);
        // }
        if (k == main_k)
        {
            // Главный процесс
            if (owner == main_k)
            {
                // Печатаем свои данные
                printed += print_array(a + block_loc * n * m, n, m, printed, max_print);
            }
            else
            {
                MPI_Status st;
                MPI_Recv(buf, n * m, MPI_DOUBLE, owner, 0 /*tag*/, comm, &st);
                printed += print_array(buf, n, m, printed, max_print);
            }
        }
        else if(owner == k)
        {
            //printf("[DEBUG][%d] Sending row number %d\n",k,block_loc);
            MPI_Send(a + block_loc * n * m, n * m, MPI_DOUBLE, main_k, 0 /*tag*/, comm);
        }
    }
}

void init_matrix(double *a, int n, int m, int p, int k,double(*f)(int i, int j, int n))
{
    /* обычные координаты, не блочные */
    int i_loc, j_loc, i_glob, j_glob, rows_block, rows;

    rows_block = get_rows_block(n, m, p, k); // кол-во блочных строк в процессе k
    //printf("[%d] Ive got %d rows\n",k,rows_block*m);
    rows = rows_block * m;

    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        
        i_glob = row_local2global(n, m, p, k, i_loc);
        //printf("[%d] For me my global %d coordinate is equal to %d global coordinate\n",k,i_loc,i_glob);
        //printf("iloc = %d, iglob = %d\n",i_loc,i_glob);
        for (j_loc = 0; j_loc < n; j_loc++)
        {
            j_glob = j_loc;
            //printf("[INIT] [%d,%d] = %lf\n",i_glob,j_glob,(*f)(i_glob, j_glob, n));
            a[i_loc * n + j_loc] = (*f)(i_glob, j_glob, n);
        }
    }
    //printf("\n\n");
}

int read_matrix(double *a, int n, int m, int p, int k, const char *name, double *buf /* n*m */)
{
    int main_k = 0;
    FILE *fp = NULL;
    int block;
    int max_blocks = (n + m - 1) / m;
    int err = 0;
    MPI_Comm comm = MPI_COMM_WORLD;

    if (k == main_k)
    {
        fp = fopen(name, "r");
        if (!fp)
            err = 1;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, comm);
    if (err > 0)
        return err;

    memset(buf, 0, n * m * sizeof(double));
    for (block = 0; block < max_blocks; block++)
    {
        int owner = block % p;
        int block_loc = block / p;

        // Главный процесс
        if (k == main_k)
        {
            err += read_array(fp, buf, n * m);

            if (owner == main_k)
            {
                // Просто скопировать буффер на место
                memcpy(a + block_loc * m * n, buf, n * m * 8);
                // for (int j = 0; j < m * n; j++)
                //     a[block_loc * m * n + j] = buf[j];
            }
            else
            {
                // Отправить буфер в процесс owner
                MPI_Send(buf, n * m, MPI_DOUBLE, owner, 0 /*tag*/, comm);
            }
        }
        else
        {
            //Остальные процессы
            if (owner == k)
            {
                MPI_Status st;
                MPI_Recv(a + block_loc * m * n, n * m, MPI_DOUBLE, main_k, 0 /*tag*/, comm, &st);
            }
        }
    }

    if (k == main_k)
    {
        fclose(fp);
        fp = 0;
    }

    MPI_Bcast(&err, 1, MPI_INT, main_k, comm);
    if (err > 0)
        return err;
    return 0;
}


int solve(double *A, double *B, int n, int m, int p, int k, double* buffer_a, double* buffer_b, double norm){
    MPI_Status st;
    int amount_of_blocks = n / m;
    int l = n % m;
    int error = 0;
    std::unique_ptr<double []> block1{new double[m*m]},block2{new double[m*m]},block3{new double[m*m]};
    std::unique_ptr<double []> a_buffer_rcv{new double[n*m]}, b_buffer_rcv{new double[n*m]};
    //Displacements for scatterv for row normalization
    std::unique_ptr<int[]> sendcounts{new int[p]};
    std::unique_ptr<int[]> displacements{new int[p]};
    double buffer = std::numeric_limits<double>::max();
    double reverse_block_norm = std::numeric_limits<double>::max();
    int reverse_block_pos = -1;
    struct
    {
        double val;
        int loc;
    } in, out;

    //Calculating displacements for row normalisation(scatterv)
    int start,stop;
    for(int i = 0; i < p; i++){
        if (i < amount_of_blocks % p) {
            start = i * (amount_of_blocks / p + 1);
            stop = start + amount_of_blocks / p;
        } else {
            start = i * (amount_of_blocks / p) + amount_of_blocks % p;
            stop = start + (amount_of_blocks / p - 1);
        }
        displacements[i] = start * m * m;
        sendcounts[i] = (stop - start + 1) * m * m;
    }
    printf("[%d] displacement = %d sendcounts = %d\n",k,displacements[k],sendcounts[k]);
    //Forward step

    for(int i = 0; i < amount_of_blocks; i++){
        buffer = std::numeric_limits<double>::max();
        reverse_block_norm = std::numeric_limits<double>::max();
        reverse_block_pos = -1;
        in.loc = -1;
        in.val = 0.0;
        //Searching for reverse blocks
        for(int j = i; j < amount_of_blocks; j++){
            if(j % p == k){
                //printf("[DEBUG] Checking block %d\n",j);
                get_block(A,n,m,j / p,i,block1.get());
                //printf("[DEBUG] Got block:\n");
                //print_matrix(block1.get(),m,m,10);
                if(inverse_matrix(block1.get(),block2.get(),m,norm) == 0){
                    //printf("[DEBUG] Got reverse block:\n");
                    //print_matrix(block2.get(),m,m,10);
                    if((buffer = norma(block2.get(),m,m)) < reverse_block_norm){
                        reverse_block_norm = buffer;
                        reverse_block_pos = j;
                    }
                }
            }
        }
        //Storing reverse norm of reverse block, so we can search with MPI_Allreduce with MPI_MAXLOC
        in.loc = reverse_block_pos;
        in.val = 1./reverse_block_norm;
        //printf("[DEBUG][%d] in.loc = %d in.val = %.3e\n",k,in.loc,in.val);
        MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
        //No reversable block found, leaving
        if(out.loc == -1){
            return -1;
        }

        //Best reversable block is misplaced, swaping it with i-row
        if(i != out.loc){
            //Rows are in different processes
            if(i % p != out.loc % p){
                
                if(k == i % p){
                    memcpy(buffer_a,A + m*n*(i / p),sizeof(double) * m * n);
                    MPI_Sendrecv(buffer_a,m*n,MPI_DOUBLE,out.loc % p,0,A + m*n*(i / p),m * n, MPI_DOUBLE, out.loc % p, 0, MPI_COMM_WORLD, &st);
                    memcpy(buffer_a,B + m*n*(i / p),sizeof(double) * m * n);
                    MPI_Sendrecv(buffer_a,m*n,MPI_DOUBLE,out.loc % p,0,B + m*n*(i / p),m * n, MPI_DOUBLE, out.loc % p, 0, MPI_COMM_WORLD, &st);
                }
                if(k == out.loc % p){
                    memcpy(buffer_a,A + m*n*(out.loc / p),sizeof(double) * m * n);
                    MPI_Sendrecv(buffer_a,m*n,MPI_DOUBLE,i % p,0,A + m*n*(out.loc / p),m * n, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);

                    memcpy(buffer_a,B + m*n*(out.loc / p),sizeof(double) * m * n);
                    MPI_Sendrecv(buffer_a,m*n,MPI_DOUBLE,i % p,0,B + m*n*(out.loc / p),m * n, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);
                }
            //Rows are in same process
            } else if (k == i % p){
                memcpy(buffer_a,A + m*n*(out.loc / p),8 * m * n);
                memcpy(A + m*n*(out.loc / p) ,A + m*n*(i / p),8 * m * n);
                memcpy(A + m*n*(i/p),buffer_a,8 * m * n);

                memcpy(buffer_a,B + m*n*(out.loc / p),8 * m * n);
                memcpy(B + m*n*(out.loc / p) ,B + m*n*(i / p),8 * m * n);
                memcpy(B + m*n*(i/p),buffer_a,8 * m * n);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //Broadcasting reverse block. We'll need it for row-normalization
        if(i % p == k){
            get_block(A,n,m,i / p,i, block1.get());
            inverse_matrix(block1.get(),block2.get(),m,norm);
        }
        MPI_Bcast(block2.get(),m*m,MPI_DOUBLE,i % p, MPI_COMM_WORLD);

        //Preparing and broadcasting the row
        if(k == i % p){
            for(int j = 0; j < amount_of_blocks; j++){
                get_block(A,n,m,i / p,j,buffer_b + m*m*j);
            }
        }
        MPI_Scatterv(buffer_b,sendcounts.get(),displacements.get(),MPI_DOUBLE,buffer_a,sendcounts[k],MPI_DOUBLE,i % p,MPI_COMM_WORLD);
        
        //Normalizing the row
        for(int j = 0; j < sendcounts[k] / (m*m); j++){
            if(displacements[k] / m / m + j < i) continue;
            multiplication(block2.get(),buffer_a + j*m*m,block3.get(),m,m,m);
            memcpy(buffer_a + j*m*m,block3.get(),m*m*8);
        }

        //Gathering row back to process i % p
        MPI_Gatherv(buffer_a,sendcounts[k],MPI_DOUBLE,buffer_b,sendcounts.get(),displacements.get(),MPI_DOUBLE,i % p, MPI_COMM_WORLD);
        if(k == i % p){
            if(l != 0){
                get_block(A,n,m,i/p,amount_of_blocks,block1.get());
                multiplication(block2.get(),block1.get(),block3.get(),m,m,l);
                put_block(A,n,m,i/p,amount_of_blocks,block3.get());
            }
            for(int j = 0; j < amount_of_blocks; j++){
                put_block(A,n,m,i/p,j,buffer_b + j*m*m);
            }
        }

        //The same operations, but for matrix B
        if(k == i % p){
            for(int j = 0; j < amount_of_blocks; j++){
                get_block(B,n,m,i / p,j,buffer_b + m*m*j);
            }
        }
        MPI_Scatterv(buffer_b,sendcounts.get(),displacements.get(),MPI_DOUBLE,buffer_a,sendcounts[k],MPI_DOUBLE,i % p,MPI_COMM_WORLD);
        //printf("[%d] Multiplication index: %d\n",k,sendcounts[k]/(m*m));
        for(int j = 0; j < sendcounts[k] / (m*m); j++){
            if(displacements[k] / m / m + j > i && displacements[k] / m / m + j != out.loc) continue;
            multiplication(block2.get(),buffer_a + j*m*m,block3.get(),m,m,m);
            memcpy(buffer_a + j*m*m,block3.get(),m*m*8);
        }
        MPI_Gatherv(buffer_a,sendcounts[k],MPI_DOUBLE,buffer_b,sendcounts.get(),displacements.get(),MPI_DOUBLE,i % p, MPI_COMM_WORLD);
        if(k == i % p){
            for(int j = 0; j < amount_of_blocks; j++){
                put_block(B,n,m,i/p,j,buffer_b + j*m*m);
            }
        }

        if(k == i % p){
            memcpy(buffer_a, A + (i / p)*m*n,m*n*8);
            memcpy(buffer_b, B + (i / p)*m*n,m*n*8);
        }

        //Broadcasting i-row
        MPI_Bcast(buffer_a,m*n,MPI_DOUBLE,i % p, MPI_COMM_WORLD);
        MPI_Bcast(buffer_b,m*n,MPI_DOUBLE,i % p, MPI_COMM_WORLD);



        //Subtracting i-row from all other rows
        for(int j = i + 1; j < amount_of_blocks; j++){
            if(j % p == k){
                get_block(A,n,m,j / p, i, block1.get());
                for(int z = i; z < amount_of_blocks; z++){
                    get_block(buffer_a,n,m,0,z,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                    get_block(A,n,m,j / p, z, block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(A,n,m,j / p, z, block2.get());
                }
                if(l != 0){
                    get_block(buffer_a,n,m,0,amount_of_blocks,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,l);
                    get_block(A,n,m,j / p, amount_of_blocks, block2.get());
                    substraction(block2.get(),block3.get(),m,l);
                    put_block(A,n,m,j / p, amount_of_blocks, block2.get());
                }
                for(int z = 0; z < i; z++){
                    get_block(buffer_b,n,m,0,z,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                    get_block(B,n,m,j / p, z, block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(B,n,m,j / p, z, block2.get());
                }
                get_block(buffer_b,n,m,0,out.loc,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                get_block(B,n,m,j / p, out.loc, block2.get());
                substraction(block2.get(),block3.get(),m,m);
                put_block(B,n,m,j / p, out.loc, block2.get());
            }
        }
        if(l != 0 && amount_of_blocks % p == k){
            get_block_raw(A,n,m,amount_of_blocks/p,i,block1.get(),l,m);
            for(int z = i; z < amount_of_blocks; z++){
                get_block(buffer_a,n,m,0,z,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),l,m,m);
                get_block_raw(A,n,m,amount_of_blocks / p, z, block2.get(),l,m);
                substraction(block2.get(),block3.get(),l,m);
                put_block_raw(A,n,m,amount_of_blocks / p, z, block2.get(),l,m);
            }
            get_block(buffer_a,n,m,0,amount_of_blocks,block2.get());
            multiplication(block1.get(),block2.get(),block3.get(),l,m,l);
            get_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block2.get(),l,l);
            substraction(block2.get(),block3.get(),l,l);
            put_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block2.get(),l,l);

            for(int z = 0; z < i; z++){
                get_block(buffer_b,n,m,0,z,block2.get());
                multiplication(block1.get(),block2.get(),block3.get(),l,m,m);
                get_block_raw(B,n,m,amount_of_blocks / p, z, block2.get(),l,m);
                substraction(block2.get(),block3.get(),l,m);
                put_block_raw(B,n,m,amount_of_blocks / p, z, block2.get(),l,m);
            }
            get_block(buffer_b,n,m,0,out.loc,block2.get());
            multiplication(block1.get(),block2.get(),block3.get(),l,m,m);
            get_block_raw(B,n,m,amount_of_blocks / p, out.loc, block2.get(),l,m);
            substraction(block2.get(),block3.get(),l,m);
            put_block_raw(B,n,m,amount_of_blocks / p, out.loc, block2.get(),l,m);
        }
    }
    //The last stage of forward step performed for bottom-border blocks
    // printf("[DEBUG] Matrix A after fowrard step:\n");
    // print_matrix(A,n,n,n);
    // printf("[DEBUG] Matrix B after fowrard step:\n");
    // print_matrix(B,n,n,n);
    if(l != 0){
        if(amount_of_blocks % p == k){ 
            get_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block1.get(),l,l);
            if(inverse_matrix(block1.get(),block2.get(),l,norm) != 0){
                error = -1;
            }
        }
        MPI_Bcast(&error,1,MPI_INT,amount_of_blocks % p, MPI_COMM_WORLD);
        if(error == -1){
            return -1;
        }
        // MPI_Bcast(block2.get(),l*l, MPI_DOUBLE, amount_of_blocks % p, MPI_COMM_WORLD);
        // if(amount_of_blocks % p == k){
        //     get_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block1.get(),l,l);
        //     multiplication(block2.get(),block1.get(),block3.get(),l,l,l);
        //     put_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block3.get(),l,l);

        //     get_block_raw(B,n,m,amount_of_blocks / p, amount_of_blocks, block1.get(),l,l);
        //     multiplication(block2.get(),block1.get(),block3.get(),l,l,l);
        //     put_block_raw(B,n,m,amount_of_blocks / p, amount_of_blocks, block3.get(),l,l);
        
        //     for(int j = 0; j < amount_of_blocks; j++){
        //         get_block_raw(B,n,m,amount_of_blocks / p, j, buffer_a + m*m*j,l,m);
        //     }
        // }
        // printf("[DEBUG] Buffer_a before scattering:\n");
        // print_matrix(buffer_a,m,amount_of_blocks * m,n); 
        // MPI_Scatterv(buffer_a, sendcounts.get(), displacements.get(), MPI_DOUBLE, buffer_b, sendcounts[k], MPI_DOUBLE, amount_of_blocks % p, MPI_COMM_WORLD);
        // printf("[DEBUG] Buffer_a after scattering:\n");
        // print_matrix(buffer_b,m,amount_of_blocks * m,n);
        // for(int j = 0; j < sendcounts[k] / (m*m); j++){
        //     get_block_raw(buffer_b + m*m*j,n,m,0,0,block1.get(),l,m);
        //     multiplication(block2.get(),block2.get(), block3.get(), l, l, m);
        //     put_block_raw(buffer_b + m*m*j,n,m,0,0,block3.get(),l,m);
        //     //memcpy(buffer_b + m*m*j,block3.get(),l*m);
        // }
        // MPI_Gatherv(buffer_b,sendcounts[k],MPI_DOUBLE,buffer_a,sendcounts.get(),displacements.get(),MPI_DOUBLE,amount_of_blocks % p, MPI_COMM_WORLD);
        // if(amount_of_blocks % p == k){
        //     for(int j = 0; j < amount_of_blocks; j++){
        //         put_block_raw(B,n,m,amount_of_blocks / p, j, buffer_a + m*m*j,l,m);
        //     }
        // }
        if(amount_of_blocks % p == k){
            get_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block1.get(),l,l);
            multiplication(block2.get(),block1.get(),block3.get(),l,l,l);
            put_block_raw(A,n,m,amount_of_blocks / p, amount_of_blocks, block3.get(),l,l);

            get_block_raw(B,n,m,amount_of_blocks / p, amount_of_blocks, block1.get(),l,l);
            multiplication(block2.get(),block1.get(),block3.get(),l,l,l);
            put_block_raw(B,n,m,amount_of_blocks / p, amount_of_blocks, block3.get(),l,l);
            for(int j = 0; j < amount_of_blocks; j++){
                get_block_raw(B,n,m,amount_of_blocks / p, j, block1.get(),l,m);
                multiplication(block2.get(),block1.get(),block3.get(),l,l,m);
                put_block_raw(B,n,m,amount_of_blocks / p, j, block3.get(),l,m);
            }
        }
    }
    
    //First stage of reverse step for bottom-border blocks
    if(l != 0){
        if(amount_of_blocks % p == k){
            memcpy(buffer_a, A + (amount_of_blocks/p)*m*n,l*n*8);
            memcpy(buffer_b, B + (amount_of_blocks/p)*m*n,l*n*8);
        }
        MPI_Bcast(buffer_a,m*n,MPI_DOUBLE,amount_of_blocks % p, MPI_COMM_WORLD);
        MPI_Bcast(buffer_b,m*n,MPI_DOUBLE,amount_of_blocks % p, MPI_COMM_WORLD);

        for(int j = 0; j < amount_of_blocks; j++){
            if(j % p == k){
                get_block(A,n,m,j / p, amount_of_blocks, block1.get());
                for(int z = 0; z < amount_of_blocks; z++){
                    get_block_raw(buffer_b,n,m,0,z,block2.get(),l,m);
                    multiplication(block1.get(),block2.get(),block3.get(),m,l,m);
                    get_block(B,n,m,j / p,z,block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(B,n,m,j / p,z,block2.get());
                }
                get_block_raw(buffer_b,n,m,0,amount_of_blocks,block2.get(),l,l);
                multiplication(block1.get(),block2.get(),block3.get(),m,l,l);
                get_block(B,n,m,j / p,amount_of_blocks,block2.get());
                substraction(block2.get(),block3.get(),m,l);
                put_block(B,n,m,j / p,amount_of_blocks,block2.get());
            }
        }
    }

    for(int i = amount_of_blocks - 1; i >= 0; i--){
        if(i % p == k){
            memcpy(buffer_a,A + (i/p)*m*n,m*n*8);
            memcpy(buffer_b,B + (i/p)*m*n,m*n*8);
        }
        MPI_Bcast(buffer_a,m*n,MPI_DOUBLE,i % p, MPI_COMM_WORLD);
        MPI_Bcast(buffer_b,n*m,MPI_DOUBLE,i % p, MPI_COMM_WORLD);
        for(int j = 0; j < i; j++){
            if(j % p == k){
                get_block(A,n,m,j/p,i,block1.get());
                for(int z = 0; z < amount_of_blocks; z++){
                    get_block(buffer_b,n,m,0,z,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,m);
                    get_block(B,n,m,j / p,z,block2.get());
                    substraction(block2.get(),block3.get(),m,m);
                    put_block(B,n,m,j / p,z,block2.get());
                }
                if(l != 0){
                    get_block(buffer_b,n,m,0,amount_of_blocks,block2.get());
                    multiplication(block1.get(),block2.get(),block3.get(),m,m,l);
                    get_block(B,n,m,j / p,amount_of_blocks,block2.get());
                    substraction(block2.get(),block3.get(),m,l);
                    put_block(B,n,m,j / p,amount_of_blocks,block2.get());
                }
            }
        }
    }
    return 0;
}
void mpi_matrix_multiplication(double *a, double *b,double *c, int n, int m, int p, int k){
    //Процессам будет необходимо сконструировать часть вектора основываясь на имеющихся у них кусках 
    //Будем считать что нужно умножить матрицы размера n*n строка-на-вектор
    //Тогда на каждом шаге процесссы соединяют куски столбцов из матрицы B в единный вектор, разсылают всем и каждый умножает строку на столбец.
    //double* peace_of_b_buffer = new double[n];
    std::unique_ptr<double[]> peace_of_b_buffer{new double[n]};
    for(int i = 0; i < n; i++){

        //Для начала, каждый процесс должен получить образец столбца на который будет домножать 
        //printf("[DEBUG][%d] Got column:\n{",k);
        for(int j = 0; j < n/m; j++){
            if(j % p == k){
                for(int z = 0; z < m; z++){
                    //printf("%lf \n",b[j/p*n*m + z*n + i]);
                    peace_of_b_buffer[j*m + z] = b[j/p*n*m + z*n + i];
                }
            }
            MPI_Bcast(&peace_of_b_buffer[j*m],m,MPI_DOUBLE,j % p,MPI_COMM_WORLD);
        }
        //printf("}\n");
        if((n % m != 0)){
            if(n/m % p == k){
                for(int z = 0; z < n % m; z++){
                    peace_of_b_buffer[n/m * m + z] = b[(n/m)/p*n*m +z*n + i];
                }
            }
            MPI_Bcast(&peace_of_b_buffer[n/m * m],n%m,MPI_DOUBLE,n/m % p,MPI_COMM_WORLD);
        }
        for(int j = 0; j < n / m; j++){
            if(j % p == k){
                for(int l = 0; l < m; l++){
                    for(int z = 0; z < n; z++){
                        c[j/p*n*m + l*n + i] += a[j/p*n*m + l*n + z] * peace_of_b_buffer[z];
                    }
                    if(j*m + l == i){
                        c[j/p*n*m + l*n + i] -= 1;
                    }
                }
            }
        }
        if((n % m != 0)&&(n/m % p == k)){
            for(int l = 0; l < n % m; l++){
                for(int z = 0; z < n; z++){
                    c[(n/m)/p*n*m + l*n + i] += a[(n/m)/p*n*m + l*n + z] * peace_of_b_buffer[z];
                }
                if((n/m)*m + l == i){
                    c[(n/m)/p*n*m + l*n + i] -= 1;
                }
            }
        }

    }
}
