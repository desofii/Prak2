#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define alpha 2.45740110027234
#define alpha2 2.47740110027234

/*void MethodProgonki (double *a, double *b, double *c, int n, double *v); // матрица размера n*n
void MethodProgonki (double *a, double *b, double *c, int n, double *v) {
    
    double *p;
    double *q;
    double *y;
    
    p = (double*)malloc(sizeof(double)*(n - 1));
    q = (double*)malloc(sizeof(double)*(n - 1));
    y = (double*)malloc(sizeof(double)*n);
    
    p[0] = - b[0] / c[0];
    q[0] = v[0] / c[0];
    
    for (int i = 0; i < n - 2; i ++) {
        p[i + 1] = - b[i + 1] / (c[i + 1] + a[i] * p[i]);
        q[i + 1] = (v[i + 1] - a[i] * q[i]) / (c[i + 1] + a[i] * p[i]);
    }
    
    y[n - 1] = (v[n - 1] - a[n - 2] * q[n - 2]) / (c[n - 1] + a[n - 2] * p[n - 2]);
    
    for (int i = n - 2; i >= 0; i --) {
        y[i] = p[i] * y[i + 1] + q[i];
    }
}*/

void MakeV (double *u1, int M, double h, double *v, double u01, double u11, double k); // функция формирует вектор правой части для вычисления следующего слоя методом прогонки
void MakeV (double *u1, int M, double h, double *v, double u01, double u11, double k) {
    
    v[0] = u01 * (M*M/2 + (M*M/2)*(-alpha*h/4 - M/(2*k) + M/2)/(alpha*h/4 - M/2 - M/(2*k))) + u1[0]*(M*M/k - M*M - M*M*M/(alpha*h - 2*M - 2*M/k) + alpha/2) + u1[1] * M * M / 2;
    v[M - 2] = M*M/2  + M*M*u11/2 + M*M*u1[M - 3]/2 + (M*M/k - M*M + alpha/2) * u1[M-2];
    for (int i = 1; i < M - 2; i ++)
        v[i] = M*M*u1[i - 1]/2 + (M*M/k - M*M + alpha/2)*u1[i] + M*M*u1[i + 1]/2;
}

void AdMethodProgonki (double *a, double *b, double *c, int n, double *v, double *y); // адаптированный под нашу задачу метод прогонки (учитываем равенство некоторых элементов нашей матрицы --- не заводим три длинных массива под элементы матрицы)
void AdMethodProgonki (double *a, double *b, double *c, int n, double *v, double *y) {
    
    double *p;
    double *q;
    
    p = (double*)malloc(sizeof(double)*(n - 1));
    q = (double*)malloc(sizeof(double)*(n - 1));
   // y = (double*)malloc(sizeof(double)*n);
    
    p[0] = - a[1] / a[0];
    q[0] = v[0] / a[0];
    
    for (int i = 0; i < n - 2; i ++) {
        p[i + 1] = - b[2] / (b[1] + b[0] * p[i]);
        q[i + 1] = (v[i + 1] - b[0] * q[i]) / (b[1] + b[0] * p[i]);
    }
    
    y[n - 1] = (v[n - 1] - c[0] * q[n - 2]) / (c[1] + c[0] * p[n - 2]);
    
    for (int i = n - 2; i >= 0; i --) {
        y[i] = p[i] * y[i + 1] + q[i]; // решение
    }
    free(p);
    free(q);
}

void Function (int M, double *a, double *b, double *c, double k, double h, double u01, double u11, double *u02, double *u12, double *u1);
void Function (int M, double *a, double *b, double *c, double k, double h, double u01, double u11, double *u02, double *u12, double *u1){
    
    double *v;
    double tmp;
    
    tmp = u1[0];
    v = (double*)malloc(sizeof(double)*(M - 1));
    MakeV(u1, M, h, v, u01, u11, k);
    AdMethodProgonki (a, b, c, M - 1, v, u1);
    *u02 = (u01*(-h*alpha/4 - M/(2*k) + M/2) - M*(u1[0] + tmp)/2) / (h*alpha/4 - M/(2*k) - M/2);
    *u12 = 1;
}

void Function1(int M, double h, int m, double *arr);
void Function1(int M, double h, int m, double *arr) {
    double *u1;
    double k;
    double u01, u11, u02, u12;
    int i;
    double a[2]; // элементы трехдиагональной матрицы
    double b[3];
    double c[2];
    
    k = 0.25;
    i = 1;
    a[0] = M*M/k + M*M - alpha/2 + 1/(4*h*h*h*(alpha*h/4 - M/2 - M/(2*k)));
    a[1] = -M*M/2;
    b[0] = -M*M/2;
    b[1] = M*M - alpha/2 + M*M/k;
    b[2] = -M*M/2;
    c[0] = -M*M/2;
    c[1] = M*M - alpha/2 + M*M/k; // задали матрицу
    
    u1 = (double*)malloc(sizeof(double)*(M - 1));
    
    u01 = 0;
    u11 = h*h*M*M*(1 - h*M);
    
    for (int i = 0; i < M - 1; i ++) {
        u1[i] = h*h*(i+1)*(i+1)*(1 - h*(i+1));
    }
    
    while(i <= 4*M*M) {
        Function(M, a, b, c, k, h, u01, u11, &u02, &u12, u1);
        u01 = u02;
        u11 = u12;
        if (i == 4*M*M*0.1) {
            printf("u0.1_0.1 = %.12f\n", u1[m - 1]);
            arr[0] = u1[m - 1];
            printf("u0.1_0.3 = %.12f\n", u1[3*m - 1]);
            arr[1] = u1[3*m - 1];
            printf("u0.1_0.5 = %.12f\n", u1[5*m - 1]);
            arr[2] = u1[5*m - 1];
            printf("u0.1_0.7 = %.12f\n", u1[7*m - 1]);
            arr[3] = u1[7*m - 1];
            printf("u0.1_0.9 = %.12f\n", u1[9*m -1]);
            arr[4] = u1[9*m - 1];
            printf("\n\n");
        }
        
        if (i == 4*M*M*0.5) {
            printf("u0.5_0.1 = %.12f\n", u1[m - 1]);
            arr[5] = u1[m - 1];
            printf("u0.5_0.3 = %.12f\n", u1[3*m - 1]);
            arr[6] = u1[3*m - 1];
            printf("u0.5_0.5 = %.12f\n", u1[5*m - 1]);
            arr[7] = u1[5*m - 1];
            printf("u0.5_0.7 = %.12f\n", u1[7*m - 1]);
            arr[8] = u1[7*m - 1];
            printf("u0.5_0.9 = %.12f\n", u1[9*m -1]);
            arr[9] = u1[9*m - 1];
            printf("\n\n");
        }
        
        if (i == 4*M*M*0.9) {
            printf("u0.9_0.1 = %.12f\n", u1[m - 1]);
            arr[10] = u1[m - 1];
            printf("u0.9_0.3 = %.12f\n", u1[3*m - 1]);
            arr[11] = u1[3*m - 1];
            printf("u0.9_0.5 = %.12f\n", u1[5*m - 1]);
            arr[12] = u1[5*m - 1];
            printf("u0.9_0.7 = %.12f\n", u1[7*m - 1]);
            arr[13] = u1[7*m - 1];
            printf("u0.9_0.9 = %.12f\n", u1[9*m -1]);
            arr[14] = u1[9*m - 1];
            printf("\n\n");
        }
            
        i ++;
    }
}

int main ()
{
    int M0, M1, M2, M3;
    int m0, m1, m2, m3;
    double h0, h1, h2, h3;
    double arr1[15], arr2[15], arr3[15];
    
    M0 = 50;
    m0 = 5;
    h0 = 0.02; // (h = 1/M)
    
    M1 = 100;
    m1 = 10;
    h1 = 0.01; // (h = 1/M)
    
    M2 = 150;
    m2 = 15;
    h2 = 0.00666666667; // (h = 1/M)
    
    M3 = 200;
    m3 = 20;
    h3 = 0.005; // (h = 1/M)
    
    Function1(M0, h0, m0, arr1);
    Function1(M1, h1, m1, arr2);
    Function1(M3, h3, m3, arr3);
    
    printf("u(50)-u(100)");
    
    for (int j = 0; j < 15; j++) {
        printf("%.12f\n", arr1[j] - arr2[j]);
    }
    
    printf("\n\n");
    printf("u(100)-u(200)");
    
    for (int j = 0; j < 15; j++) {
        printf("%.12f\n", arr2[j] - arr3[j]);
    }
    
    printf("\n\n");
    printf("u(50)-u(100) / u(100)-u(200)");
    
    for (int j = 0; j < 15; j++) {
        printf("%.12f\n", (arr1[j] - arr2[j]) / (arr2[j] - arr3[j]));
    }
    return 0;
}
