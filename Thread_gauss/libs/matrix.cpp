#include "matrix.h"
#include <stdio.h>
#include <math.h>
void print_matrix(double* A, int n, int m, int r){
    int i,j;
    int nn = (n < r)? n : r;
    int mm = (m < r)? m : r; 
    for(i = 0; i < nn; i++){
        for(j = 0; j < mm; j++){
            printf("%3e ",A[i*m + j]);
        }
        printf("\n");
    }
    printf("\n");
}
double norma(double *a, int p, int q)  // l1
{
  int i = 0, j = 0;
  double sum = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
  // if (p % 3 == 0 && q % 3 == 0) {
  //   for (i = 0; i < p; i++) {
  //     for (j = 0; j < q; j += 3) {
  //       sum1 += (a[i * q + j] * a[i * q + j]) +
  //               (a[i * q + j + 1] * a[i * q + j + 1]) +
  //               (a[i * q + j + 1] * a[i * q + j + 1]);

  //       i++;
  //       sum2 += (a[i * q + j] * a[i * q + j]) +
  //               (a[i * q + j + 1] * a[i * q + j + 1]) +
  //               (a[i * q + j + 1] * a[i * q + j + 1]);
  //       i++;
  //       sum3 += (a[i * q + j] * a[i * q + j]) +
  //               (a[i * q + j + 1] * a[i * q + j + 1]) +
  //               (a[i * q + j + 1] * a[i * q + j + 1]);
  //     }
  //     sum = sum1 + sum2 + sum3;
  //     sum1 = 0.0;
  //     sum2 = 0.0;
  //     sum3 = 0.0;
  //   }
  // } else {
    for (i = 0; i < p; i++) {
      for (j = 0; j < q; j++) {
        sum += (a[i * q + j]) * (a[i * q + j]);
      }
    }
  //}

  return sum;
}
int inverse_matrix(double *a, double *x, int n, double norm) {
  int i = 0, j = 0, k = 0;
  double eps = 1e-15;
  double tmp1,tmp2,tmp3;
  double leader_norm;
  int leader_pos;
  if (n % 3 == 0) {
    for (i = 0; i < n; i += 3)
      for (j = 0; j < n; j += 3) {
        if (i == j) {
          x[i * n + i] = 1.0;
          x[i * n + j + 1] = 0.0;
          x[i * n + j + 2] = 0.0;

          x[(i + 1) * n + j] = 0.0;
          x[(i + 1) * n + j + 1] = 1.0;

          x[(i + 1) * n + j + 2] = 0.0;

          x[(i + 2) * n + j] = 0.0;
          x[(i + 2) * n + j + 1] = 0.0;
          x[(i + 2) * n + j + 2] = 1.0;

        } else {
          x[i * n + j] = 0.0;
          x[i * n + j + 1] = 0.0;
          x[i * n + j + 2] = 0.0;

          x[(i + 1) * n + j] = 0.0;
          x[(i + 1) * n + j + 1] = 0.0;
          x[(i + 1) * n + j + 2] = 0.0;

          x[(i + 2) * n + j] = 0.0;
          x[(i + 2) * n + j + 1] = 0.0;
          x[(i + 2) * n + j + 2] = 0.0;
        }
      }
  } else {
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++) x[i * n + j] = (i == j);
  }

  if (n % 3 == 0) {
    for (k = 0; k < n; k++) {
      leader_norm = fabs(a[k*n + k]);
      leader_pos = k;
      for(i = k; i < n; i++){
        if(fabs(a[i*n + k]) > leader_norm){
          leader_norm = fabs(a[i*n + k]);
          leader_pos = i;
        }
      }
      for(i = 0; i < n; i += 3){
          tmp1 = a[k*n + i];
          tmp2 = a[k*n + i + 1];
          tmp3 = a[k*n + i + 2];
          a[k*n + i] = a[leader_pos*n + i];
          a[k*n + i + 1] = a[leader_pos*n + i + 1];
          a[k*n + i + 2] = a[leader_pos*n + i + 2];
          a[leader_pos*n + i] = tmp1;
          a[leader_pos*n + i + 1] = tmp2;
          a[leader_pos*n + i + 2] = tmp3;

          tmp1 = x[k*n + i];
          tmp2 = x[k*n + i + 1];
          tmp3 = x[k*n + i + 2];
          x[k*n + i] = x[leader_pos*n + i];
          x[k*n + i + 1] = x[leader_pos*n + i + 1];
          x[k*n + i + 2] = x[leader_pos*n + i + 2];
          x[leader_pos*n + i] = tmp1;
          x[leader_pos*n + i + 1] = tmp2;
          x[leader_pos*n + i + 2] = tmp3;
      }
      if (fabs(a[k * n + k]) <= eps * norm)  // if a_ss == 0
      {
        return -1;
      }
      for (i = 0; i < n; i += 3) {
        x[k * n + i] /= a[k * n + k];
        x[k * n + i + 1] /= a[k * n + k];
        x[k * n + i + 2] /= a[k * n + k];
      }

      for (i = 0; i < n - k; i++) a[k * n + (n - i - 1)] /= a[k * n + k];

      for (i = 0; i < n; i++) {
        if (i == k) continue;
        for (j = 0; j < n; j += 3) {
          x[i * n + j] -= x[k * n + j] * a[i * n + k];
          x[i * n + j + 1] -= x[k * n + j + 1] * a[i * n + k];
          x[i * n + j + 2] -= x[k * n + j + 2] * a[i * n + k];
        }
        for (j = 0; j < n - k; j++) {
          a[i * n + (n - j - 1)] -= a[k * n + (n - j - 1)] * a[i * n + k];
        }
      }
    }

  } else {
    for (k = 0; k < n; k++) {
      // for (i = k; i < n; i++) {
      //   if (fabs(a[k * n + k]) <= eps * norm)  // if a_ss == 0
      //   {
      //     return -1;
      //   }
      // }
      leader_pos = k;
      leader_norm = fabs(a[k*n + k]);
      for(i = k; i < n; i++){
        if(fabs(a[i*n + k]) > leader_norm){
          leader_norm = fabs(a[i*n + k]);
          leader_pos = i;
        }
      }
      for(i = 0; i < n; i++){
        tmp1 = a[k*n + i];
        a[k*n + i] = a[leader_pos*n + i];
        a[leader_pos*n + i] = tmp1;

        tmp1 = x[k*n + i];
        x[k*n + i] = x[leader_pos*n + i];
        x[leader_pos*n + i] = tmp1;
      }
      if(fabs(a[k*n + k]) <= eps * norm){
        return -1;
      }
      for (i = 0; i < n; i++) x[k * n + i] /= a[k * n + k];

      for (i = 0; i < n - k; i++) a[k * n + (n - i - 1)] /= a[k * n + k];

      for (i = 0; i < n; i++) {
        if (i == k) continue;
        for (j = 0; j < n; j++) x[i * n + j] -= x[k * n + j] * a[i * n + k];
        for (j = 0; j < n - k; j++)
          a[i * n + (n - j - 1)] -= a[k * n + (n - j - 1)] * a[i * n + k];
      }
    }
  }

  return 0;
}
void put_block(double* A, int n, int m, int i, int j, double* block){

    int p,q,t,s;
    if(i != (n / m)){
        p = m;
    } else {
        p = n % m;
    }
    if(j != (n/m)){
        q = m;
    } else {
        q = n % m;
    }
    // p,q - block sizes
    if(p % 3 == 0 && q % 3 == 0){
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t += 3){
                A[i*m*n + j*m + s*n + t] = block[s * q + t];
                A[i*m*n + j*m + s*n + t + 1] = block[s * q + t + 1];
                A[i*m*n + j*m + s*n + t + 2] = block[s * q + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                A[i*m*n + j*m + s*n + t] = block[s * q + t];
                A[i*m*n + j*m + s*n + t + 1] = block[s * q + t + 1];
                A[i*m*n + j*m + s*n + t + 2] = block[s * q + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                A[i*m*n + j*m + s*n + t] = block[s * q + t];
                A[i*m*n + j*m + s*n + t + 1] = block[s * q + t + 1];
                A[i*m*n + j*m + s*n + t + 2] = block[s * q + t + 2];
            }
        }
    } else {
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t++){
                A[i*m*n + j*m + s*n + t] = block[s*q + t];
            }
        }
    }
    // printf("[DEBUG] Got block:\n");
    // print_matrix(block,p,q);

}
void get_block(double* A, int n, int m, int i, int j, double* block){

    int p,q,t,s;
    if(i != (n / m)){
        p = m;
    } else {
        p = n % m;
    }
    if(j != (n/m)){
        q = m;
    } else {
        q = n % m;
    }
    // p,q - block sizes
    if(p % 3 == 0 && q % 3 == 0){
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t += 3){
                block[s * q + t] = A[i*m*n + j*m + s*n + t];
                block[s * q + t + 1] = A[i*m*n + j*m + s*n + t + 1];
                block[s * q + t + 2] = A[i*m*n + j*m + s*n + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                block[s * q + t] = A[i*m*n + j*m + s*n + t];
                block[s * q + t + 1] = A[i*m*n + j*m + s*n + t + 1];
                block[s * q + t + 2] = A[i*m*n + j*m + s*n + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                block[s * q + t] = A[i*m*n + j*m + s*n + t];
                block[s * q + t + 1] = A[i*m*n + j*m + s*n + t + 1];
                block[s * q + t + 2] = A[i*m*n + j*m + s*n + t + 2];
            }
        }
    } else {
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t++){
                block[s*q + t] = A[i*m*n + j*m + s*n + t];
            }
        }
    }
    // printf("[DEBUG] Got block:\n");
    // print_matrix(block,p,q);

}

void put_block_raw(double* A, int n, int m, int i, int j, double* block, int p, int q){

    int /*p,q,*/t,s;
    // if(i != (n / m)){
    //     p = m;
    // } else {
    //     p = n % m;
    // }
    // if(j != (n/m)){
    //     q = m;
    // } else {
    //     q = n % m;
    // }
    // p,q - block sizes
    if(p % 3 == 0 && q % 3 == 0){
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t += 3){
                A[i*m*n + j*m + s*n + t] = block[s * q + t];
                A[i*m*n + j*m + s*n + t + 1] = block[s * q + t + 1];
                A[i*m*n + j*m + s*n + t + 2] = block[s * q + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                A[i*m*n + j*m + s*n + t] = block[s * q + t];
                A[i*m*n + j*m + s*n + t + 1] = block[s * q + t + 1];
                A[i*m*n + j*m + s*n + t + 2] = block[s * q + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                A[i*m*n + j*m + s*n + t] = block[s * q + t];
                A[i*m*n + j*m + s*n + t + 1] = block[s * q + t + 1];
                A[i*m*n + j*m + s*n + t + 2] = block[s * q + t + 2];
            }
        }
    } else {
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t++){
                A[i*m*n + j*m + s*n + t] = block[s*q + t];
            }
        }
    }
    // printf("[DEBUG] Got block:\n");
    // print_matrix(block,p,q);

}
void get_block_raw(double* A, int n, int m, int i, int j, double* block, int p, int q){

    int /*p,q,*/t,s;
    // if(i != (n / m)){
    //     p = m;
    // } else {
    //     p = n % m;
    // }
    // if(j != (n/m)){
    //     q = m;
    // } else {
    //     q = n % m;
    // }
    // p,q - block sizes
    if(p % 3 == 0 && q % 3 == 0){
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t += 3){
                block[s * q + t] = A[i*m*n + j*m + s*n + t];
                block[s * q + t + 1] = A[i*m*n + j*m + s*n + t + 1];
                block[s * q + t + 2] = A[i*m*n + j*m + s*n + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                block[s * q + t] = A[i*m*n + j*m + s*n + t];
                block[s * q + t + 1] = A[i*m*n + j*m + s*n + t + 1];
                block[s * q + t + 2] = A[i*m*n + j*m + s*n + t + 2];
            }
            s++;
            for(t = 0; t < q; t += 3){
                block[s * q + t] = A[i*m*n + j*m + s*n + t];
                block[s * q + t + 1] = A[i*m*n + j*m + s*n + t + 1];
                block[s * q + t + 2] = A[i*m*n + j*m + s*n + t + 2];
            }
        }
    } else {
        for(s = 0; s < p; s++){
            for(t = 0; t < q; t++){
                block[s*q + t] = A[i*m*n + j*m + s*n + t];
            }
        }
    }
    // printf("[DEBUG] Got block:\n");
    // print_matrix(block,p,q);

}

double formula(int n, int x, int y, int mode){
    switch (mode)
    {
    case 1:
        if(x > y){
            return (n - x);
        } else {
            return (n - y);
        }
        break;
    case 2:
        if(x > y){
            return x + 1.0;
        } else {
            return y + 1.0;
        }
        break;
    case 3:
        if(x > y){
            return (x - y);
        } else {
            return (y - x);
        }
        break;
    case 4:
        return 1.0 / (x + y + 1.0);
        break;
    }
    return -1;
}


void ones(double *a, int n){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      a[i*n + j] = (i == j);
    }
  }
}

void multiplication(double *a, double *b, double *c, int n, int m,
                    int q)  // n*m m*q = n*q
{
  int i = 0, j = 0, k = 0;
  double x11, x12, x13, x21, x22, x23, x31, x32, x33;

  if (n % 3 == 0 && m % 3 == 0 && q % 3 == 0) {
    for (i = 0; i < n; i += 3)
      for (j = 0; j < q; j += 3) {
        x11 = 0;
        x12 = 0;
        x13 = 0;

        x21 = 0;
        x22 = 0;
        x23 = 0;

        x31 = 0;
        x32 = 0;
        x33 = 0;
        c[i * q + j] = 0;
        c[i * q + j + 1] = 0;
        c[i * q + j + 2] = 0;

        c[(i + 1) * q + j] = 0;
        c[(i + 1) * q + j + 1] = 0;
        c[(i + 1) * q + j + 2] = 0;

        c[(i + 2) * q + j] = 0;
        c[(i + 2) * q + j + 1] = 0;
        c[(i + 2) * q + j + 2] = 0;
        // c[i * q + j] = 0;
        for (k = 0; k < m; k += 3) {
          x11 += a[i * m + k] * b[k * q + j];
          x12 += a[i * m + k] * b[k * q + j + 1];
          x13 += a[i * m + k] * b[k * q + j + 2];
          x21 += a[(i + 1) * m + k] * b[k * q + j];
          x22 += a[(i + 1) * m + k] * b[k * q + j + 1];
          x23 += a[(i + 1) * m + k] * b[k * q + j + 2];
          x31 += a[(i + 2) * m + k] * b[k * q + j];
          x32 += a[(i + 2) * m + k] * b[k * q + j + 1];
          x33 += a[(i + 2) * m + k] * b[k * q + j + 2];

          x11 += a[i * m + k + 1] * b[(k + 1) * q + j];
          x12 += a[i * m + k + 1] * b[(k + 1) * q + j + 1];
          x13 += a[i * m + k + 1] * b[(k + 1) * q + j + 2];
          x21 += a[(i + 1) * m + k + 1] * b[(k + 1) * q + j];
          x22 += a[(i + 1) * m + k + 1] * b[(k + 1) * q + j + 1];
          x23 += a[(i + 1) * m + k + 1] * b[(k + 1) * q + j + 2];
          x31 += a[(i + 2) * m + k + 1] * b[(k + 1) * q + j];
          x32 += a[(i + 2) * m + k + 1] * b[(k + 1) * q + j + 1];
          x33 += a[(i + 2) * m + k + 1] * b[(k + 1) * q + j + 2];

          x11 += a[i * m + k + 2] * b[(k + 2) * q + j];
          x12 += a[i * m + k + 2] * b[(k + 2) * q + j + 1];
          x13 += a[i * m + k + 2] * b[(k + 2) * q + j + 2];
          x21 += a[(i + 1) * m + k + 2] * b[(k + 2) * q + j];
          x22 += a[(i + 1) * m + k + 2] * b[(k + 2) * q + j + 1];
          x23 += a[(i + 1) * m + k + 2] * b[(k + 2) * q + j + 2];
          x31 += a[(i + 2) * m + k + 2] * b[(k + 2) * q + j];
          x32 += a[(i + 2) * m + k + 2] * b[(k + 2) * q + j + 1];
          x33 += a[(i + 2) * m + k + 2] * b[(k + 2) * q + j + 2];
        }
        c[i * q + j] += x11;
        c[i * q + j + 1] += x12;
        c[i * q + j + 2] += x13;
        c[(i + 1) * q + j] += x21;
        c[(i + 1) * q + j + 1] += x22;
        c[(i + 1) * q + j + 2] += x23;
        c[(i + 2) * q + j] += x31;
        c[(i + 2) * q + j + 1] += x32;
        c[(i + 2) * q + j + 2] += x33;
      }
  } else {
    for (i = 0; i < n; i++)
      for (j = 0; j < q; j++) {
        c[i * q + j] = 0;
        for (k = 0; k < m; k++) {
          c[i * q + j] += a[i * m + k] * b[k * q + j];
        }
      }
  }
}
void substraction(double *a, double *b, int n, int m)  // -
{
  int i = 0, j = 0, k = 0, t = 0;
  if (n % 3 == 0 && m % 3 == 0) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < m; j += 3) {
        a[i * m + j] = a[i * m + j] - b[i * m + j];
        a[i * m + j + 1] = a[i * m + j + 1] - b[i * m + j + 1];
        a[i * m + j + 2] = a[i * m + j + 2] - b[i * m + j + 2];
      }
      i++;
      for (k = 0; k < m; k += 3) {
        a[i * m + k] = a[i * m + k] - b[i * m + k];
        a[i * m + k + 1] = a[i * m + k + 1] - b[i * m + k + 1];
        a[i * m + k + 2] = a[i * m + k + 2] - b[i * m + k + 2];
      }
      i++;
      for (t = 0; t < m; t += 3) {
        a[i * m + t] = a[i * m + t] - b[i * m + t];
        a[i * m + t + 1] = a[i * m + t + 1] - b[i * m + t + 1];
        a[i * m + t + 2] = a[i * m + t + 2] - b[i * m + t + 2];
      }
    }

  } else {
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++) a[i * m + j] = a[i * m + j] - b[i * m + j];
  }
}

void zero(double *a, int n, int m){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      a[i*m + j] = 0.0;
    }
  }
}

int read_matrix(int n, int s, double *a, double *b, FILE *fin) {
    int i = 0, j = 0, p;

    if (s == 0) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                if (fscanf(fin, "%lf", &a[i * n + j]) != 1)
                    return -1;
    } else {
        if (s == 1) {
            if (n % 3 == 0) {

                for (i = 0; i < n; i += 3) {
                    for (j = 0; j < n; j += 3) {
                        if (i > j) {
                            a[i * n + j] = double(n - i);
                            a[(i + 1) * n + j + 1] = double(n - i - 1);
                            a[(i + 2) * n + j + 2] = double(n - i - 2);
                        } else {
                            a[i * n + j] = double(n - j);
                            a[(i + 1) * n + j + 1] = double(n - j - 1);
                            a[(i + 2) * n + j + 2] = double(n - j - 2);
                        }

                        if (i > j + 1)
                            a[i * n + j + 1] = double(n - i);
                        else a[i * n + j + 1] = double(n - j - 1);
                        if (i > j + 2)
                            a[i * n + j + 2] = double(n - i);
                        else a[i * n + j + 2] = double(n - j - 2);

                        if (i + 1 > j)
                            a[(i + 1) * n + j] = double(n - i - 1);
                        else a[(i + 1) * n + j] = double(n - j);

                        if (i + 1 > j + 2)
                            a[(i + 1) * n + j + 2] = double(n - i - 1);
                        else a[(i + 1) * n + j + 2] = double(n - j - 2);


                        if (i + 2 > j)
                            a[(i + 2) * n + j] = double(n - i - 2);
                        else a[(i + 2) * n + j] = double(n - j);
                        if (i + 2 > j + 1)
                            a[(i + 2) * n + j + 1] = double(n - i - 2);
                        else a[(i + 2) * n + j + 1] = double(n - j - 1);


                    }
                }
            } else {
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++) {
                        if (i > j)
                            a[i * n + j] = double(n - i);
                        else a[i * n + j] = double(n - j);

                    }
            }
        } else if (s == 2) {
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++) {
                    if (i > j)
                        a[i * n + j] = double(i + 1.0);
                    else a[i * n + j] = double(j + 1.0);

                }
        } else if (s == 3) {
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++) {
                    if (i > j)
                        a[i * n + j] = double(i - j);
                    else a[i * n + j] = double(j - i);
                }
        } else if (s == 4) {
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++) {
                    a[i * n + j] = double(1.0 / (i + j + 1.0));

                }
        } else return -2;
    }

    for (i = 0; i < n; i++) {
        b[i] = 0;
        for (j = 0; j < (n + 1) / 2; j++) {

            b[i] += a[i * n + 2 * j];
        }


    }

    return 0;
}
