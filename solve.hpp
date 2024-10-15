#pragma once
#include <iostream>
#include <math.h>
#include <float.h>
#include "mult.hpp"
#include <iomanip>
#include <cstring>
//#include "reader.hpp"
#define UNUSED(x) (void)(x)
#define eps -DBL_MAX

using namespace std;

void get_block (double *matr, double *block , int n, int m, int i , int j )
{  
    for (int i=0; i < m*m; i++) 
        block[i] = 0;
    int ii=0,jj=0;
    if(m == 0 ) {m = 1;}
    int k = n / m;
    int l = n - k * m;
    int v = (j < k ? m : l), h = (i < k ? m : l);
    double *ma_bl = matr + i * n * m+ j * m;
    for (ii = 0; ii < h; ii++)
    {
        for (jj = 0; jj < v; jj++)
        {
            block[ii *m+jj]=ma_bl[ii *n+jj];
        }   
    }
}
void put_block (double *matr, double *block , int n, int m, int i , int j )
{
int ii=0,jj=0;
int k = n / m;
int l = n - k * m;
int v = (j < k ? m : l), h = (i < k ? m : l);
double *ma_bl = matr + i * n * m+ j * m;
for (ii = 0; ii < h; ii++)
{
    for (jj = 0; jj < v; jj++)
    {
        ma_bl[ii *n+jj]=block[ii *m+jj];
    }   
}
}

double norma(double* block, int m) {
    int i,j;
    double norm_m=0, temp;
    for(i=0;i<m;i++)  
    {
        temp=0;
        for(j=0;j<m;j++)
            temp+= fabs(block[j*m + i]);
        if(temp>norm_m)
            norm_m=temp;
    }
    return norm_m;
} 

int treug(double * a, double * b, int n, double norma, double* c) {
    int i;
    int j;
    int k;
    int t;
    double p;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i * n + j] = 0;
        }
        for (int i = 0; i < n; i++) {
            b[i * n + i] = 1;
        }
    }

    for (i = 0; i < n; i++) {
        t = -1;
        for (j = i; j < n; j++)
            if (a[j * n + i] >  5e-15 * norma || a[j * n + i] < -5e-15 * norma) {
                //printf("%e\n", a[j * n + i]);
                t = 1;
            }
        if (t == -1) {
            return -1;
        }

        p = a[i * n + i];
        t = i;

        for (j = i; j < n; j++) {
            if (fabs(a[j * n + i]) > fabs(p)) {
                p = a[j * n + i];
                t = j;
            }
        }

        for (j = 0; j < n; j++) c[j] = a[t * n + j];
        for (j = 0; j < n; j++) a[t * n + j] = a[i * n + j];
        for (j = 0; j < n; j++) a[i * n + j] = c[j];

        for (j = 0; j < n; j++) c[j] = b[t * n + j];
        for (j = 0; j < n; j++) b[t * n + j] = b[i * n + j];
        for (j = 0; j < n; j++) b[i * n + j] = c[j];

        p = a[i * n + i];

        for (j = 0; j < n; j++) {
            a[i * n + j] = a[i * n + j] / p;
            b[i * n + j] = b[i * n + j] / p;
        }

        for (j = i + 1; j < n; j++) {
            p = a[j * n + i];
            for (k = 0; k < n; k++) {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }

    }
    return 1;
}

void diag(double * a, double * b, int n) {
    int i, j, k;
    double p;
    for (i = n - 1; i > 0; i--) {
        for (j = 0; j < i; j++) {
            p = a[j * n + i];
            for (k = 0; k < n; k++) {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }
    }
}

double ravno(double x, double norma)
{
    if(fabs(x)<1e-7*norma){
        return 0.0;}
    else {return x;}
}

void ravnoBlock(double *block, int m, double norma)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            block[j*m + i] = ravno(block[j*m + i], norma);
        }
    }
}
void pcord(int i, int j)
{
    printf("(%d, %d)\n", i, j);
}

int findmax(double* matr, double* block, int n, int m,int l, int j, double norm, double* buffer) {
    int k = -1;
    int t;
    double tmp;
    double max_n = DBL_MAX;
    double* inverse = new double[m * m];
    for (int i = l; i < n/m; i++) {
        get_block(matr, block, n ,m , i, j);
        t = treug(block, inverse, m, norm, buffer);
        if (t==-1) {
            continue;
        }
        diag(block,inverse,m);
        tmp = norma(inverse, m);
        if (tmp < max_n) {
            max_n = tmp;
            k = i;
        }
        
    }
    delete [] inverse;
    return k;
}


void swap_rows(double* matr, int k, int l, int n, int m) {
    double* tmp = new double[m*n];
    for (int i = 0; i < m*n; i++) {
        tmp[i] = matr[k*m * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*k * n + i] = matr[m*l * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*l * n + i] = tmp[i];
    }
    delete[] tmp;
}

void subtraction(double* a, double* b, int m) {
    for (int i = 0; i < m*m; i++) {
        a[i] -= b[i];
    }    
}
void initialize_matrix(double* solution, int n) {
    std::memset(solution, 0, n * n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        solution[i * n + i] = 1.0;
    }
}
void subtraction(double *a, double *b, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            a[i * col + j] -= b[i * col + j];
        }
    }
}

bool inverse_matrixx(double *matrix, double *inverse_matrix, int n, int m, double norm)
{
    int i, j, k, t;
    double p, *a, *b, **a_rows, **b_rows, *temp_row;

    a = new double[n * n];
    b = new double[n * n];
    a_rows = new double*[n];
    b_rows = new double*[n];

    for (i = 0; i < n; i++)
    {
        a_rows[i] = &a[i * n];
        b_rows[i] = &b[i * n];
        for (j = 0; j < n; j++)
        {
            a_rows[i][j] = matrix[i * m + j];
            b_rows[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (i = 0; i < n; i++)
    {
        t = -1;
        p = 0.0;
        for (j = i; j < n; j++)
        {
            double val = fabs(a_rows[j][i]);
            if (val > fabs(p))
            {
                p = a_rows[j][i];
                t = j;
            }
        }

        if (t == -1 || fabs(p) <= 5e-15 * norm)
        {
            delete[] a;
            delete[] b;
            delete[] a_rows;
            delete[] b_rows;
            return false;
        }

        if (t != i)
        {
            temp_row = a_rows[t];
            a_rows[t] = a_rows[i];
            a_rows[i] = temp_row;

            temp_row = b_rows[t];
            b_rows[t] = b_rows[i];
            b_rows[i] = temp_row;
        }

        double p_inv = 1.0 / a_rows[i][i];
        a_rows[i][i] = 1.0;
        for (j = i + 1; j < n; j++)
        {
            a_rows[i][j] *= p_inv;
        }
        for (j = 0; j < n; j++)
        {
            b_rows[i][j] *= p_inv;
        }

        for (j = i + 1; j < n; j++)
        {
            p = a_rows[j][i];
            a_rows[j][i] = 0.0;
            for (k = i + 1; k < n; k++)
            {
                a_rows[j][k] -= a_rows[i][k] * p;
            }
            for (k = 0; k < n; k++)
            {
                b_rows[j][k] -= b_rows[i][k] * p;
            }
        }
    }

    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < i; j++)
        {
            p = a_rows[j][i];
            a_rows[j][i] = 0.0;
            for (k = 0; k < n; k++)
            {
                b_rows[j][k] -= b_rows[i][k] * p;
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            inverse_matrix[i * m + j] = b_rows[i][j];
        }
    }

    delete[] a;
    delete[] b;
    delete[] a_rows;
    delete[] b_rows;
    return true;
}

void change_block_place(double *matrix, int n, int m, int from_i, int from_j, int to_i, int to_j)
{
    // меняем строки
    int start_from = from_i * m;
    int start_to = to_i * m;
    double tmp;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tmp = matrix[(start_from + i) * n + j];
            matrix[(start_from + i) * n + j] = matrix[(start_to + i) * n + j];
            matrix[(start_to + i) * n + j] = tmp;
        }
    }

    // меняем столбцы
    start_from = from_j * m;
    start_to = to_j * m;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            tmp = matrix[i * n + (start_from + j)];
            matrix[i * n + (start_from + j)] = matrix[i * n + (start_to + j)];
            matrix[i * n + (start_to + j)] = tmp;
        }
    }
}

double matrix_norm(double *matrix, int n, int m)
{
    double norm = 0;
    for (int j = 0; j < m; j++)
    {
        double buf = 0;
        for (int i = 0; i < n; i++)
        {
            buf += fabs(matrix[i * n + j]);
        }
        if (norm < buf)
        { // возможно заменить на norm - buf < 0
            norm = buf;
        }
    }
    return norm;
}
int solve(double *matrix, double *inverse_matrix, double *block,
                                    double *inv_block, double *help_block, int *permutation,
                                    int n, int m, double matrix_norm)
{
    int block_count = n / m;                 // количество блоков m на m
    int remainder = n - block_count * m;     // размер для прямоугольных блоков
    int block_limit = (remainder != 0 ? block_count + 1 : block_count); // предел индексирования
    int min_index_i;                         // индекс блока с минимальной обратной
    int min_index_j;                         // индекс блока с минимальной обратной

    // инициализируем строку permutation
    for (int i = 0; i < block_count; i++)
    {
        permutation[i] = i;
    }

    for (int step = 0; step < block_limit; step++) // step - номер шага алгоритма
    {
        double min_inverse_norm = __DBL_MAX__; // минимальная норма обратного блока
        int inverse_block_count = 0;           // счетчик для блоков с обратной матрицей
        min_index_i = step;
        min_index_j = step;

        // Находим блок с минимальной нормой обратной и сдвигаем его в верхний левый угол подматрицы
        for (int i = step; i < block_count; i++) // (i , j) - индекс блока
        {
            get_block(matrix, block, n, m, i, step);
            if (inverse_matrixx(block, block, m, m, matrix_norm))
            {
                inverse_block_count++;
                double buffer = norma(block, m);
                if (buffer < min_inverse_norm)
                {
                    min_inverse_norm = buffer;
                    min_index_i = i;
                    for (int ii = 0; ii < m; ii++)
                    {
                        for (int jj = 0; jj < m; jj++)
                        {
                            inv_block[ii * m + jj] = block[ii * m + jj]; // обратный блок с минимальной нормой
                        }
                    }
                }
            }
        }

        if (step == block_count)
        {
            get_block(matrix, block, n, m, step, step);
            if (inverse_matrixx(block, block, remainder, m, matrix_norm))
            {
                inverse_block_count++;
                min_index_i = step;
                min_index_j = step;
                for (int ii = 0; ii < remainder; ii++)
                {
                    for (int jj = 0; jj < remainder; jj++)
                    {
                        inv_block[ii * m + jj] = block[ii * m + jj];
                    }
                }
            }
        }

        if (inverse_block_count == 0)
        {
            printf("ERROR in Solve Function: No block with inverse!\n");
            return 1;
        }
        if (step != block_count)
        {
            change_block_place(matrix, n, m, min_index_i, min_index_j, step, step);
            change_block_place(inverse_matrix, n, m, min_index_i, step, step, step);

            int temp = permutation[step];
            permutation[step] = permutation[min_index_j];
            permutation[min_index_j] = temp;
        }

        // домножаем на обратную в нужной строке и потом вычитаем
        for (int i = step; i < block_limit; i++)
        {
            for (int j = step; j < block_limit; j++)
            {
                if (i == step)
                {
                    get_block(matrix, block, n, m, i, j);
                    if (j == step)
                    {
                        for (int ii = 0; ii < m; ii++)
                        {
                            for (int jj = 0; jj < m; jj++)
                            {
                                block[ii * m + jj] = (ii == jj ? 1 : 0);
                            }
                        }

                        put_block(matrix, block, n, m, i, j);

                        for (int jj = 0; jj < block_limit; jj++)
                        {
                            get_block(inverse_matrix, block, n, m, i, jj);
                            mult(inv_block, block, help_block, m, m, m, m, matrix_norm);
                            put_block(inverse_matrix, help_block, n, m, i, jj);
                        }

                        for (int jj = step + 1; jj < block_limit; jj++)
                        {
                            get_block(matrix, block, n, m, i, jj);
                            mult(inv_block, block, help_block, m, m, m, m, matrix_norm);
                            put_block(matrix, help_block, n, m, i, jj);
                        }

                        continue;
                    }
                }
                else
                {
                    if (j == step)
                    {
                        get_block(matrix, inv_block, n, m, i, step);
                        for (int jj = 0; jj < block_limit; jj++)
                        {
                            get_block(inverse_matrix, help_block, n, m, step, jj);
                            mult(inv_block, help_block, block, m, m, m, m, matrix_norm);
                            get_block(inverse_matrix, help_block, n, m, i, jj);
                            subtraction(help_block, block, m, m);
                            put_block(inverse_matrix, help_block, n, m, i, jj);
                        }
                    }

                    get_block(matrix, help_block, n, m, step, j);
                    mult(inv_block, help_block, block, m, m, m, m, matrix_norm);
                    get_block(matrix, help_block, n, m, i, j);
                    subtraction(help_block, block, m, m);
                    put_block(matrix, help_block, n, m, i, j);
                }
            }
        }
    }

    // обратный ход
    for (int i = block_limit - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            get_block(matrix, inv_block, n, m, j, i);

            for (int jj = 0; jj < block_limit; jj++)
            {
                get_block(inverse_matrix, help_block, n, m, i, jj);
                mult(inv_block, help_block, block, m, m, m, m, matrix_norm);
                get_block(inverse_matrix, help_block, n, m, j, jj);
                subtraction(help_block, block, m, m);
                put_block(inverse_matrix, help_block, n, m, j, jj);
            }
        }
    }
    return 0;
}
