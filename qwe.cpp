#pragma once
#include <cmath>
#include <iostream>
#include <limits.h>
#include <stdio.h>

#include "LN-matrix-output.h"

void get_block(double *matrix, double *block, int n, int m, int i, int j)
{
    int k = n / m;
    int l = n - k * m;
    int v = (i < k ? m : l);
    int h = (j < k ? m : l);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < m; jj++)
        {
            block[ii * m + jj] = 0;
        }
    }

    double *matrix_block = matrix + i * n * m + j * m;
    for (int ii = 0; ii < v; ii++)
    {
        for (int jj = 0; jj < h; jj++)
        {
            block[ii * m + jj] = matrix_block[ii * n + jj];
        }
    }
}

void put_block(double *matrix, double *block, int n, int m, int i, int j)
{
    int k = n / m;
    int l = n - k * m;
    int v = (j < k ? m : l);
    int h = (i < k ? m : l);
    double *matrix_block = matrix + i * n * m + j * m;
    for (int ii = 0; ii < h; ii++)
    {
        for (int jj = 0; jj < v; jj++)
        {
            matrix_block[ii * n + jj] = block[ii * m + jj];
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

void matrix_mult(double *a, double *b, double *c, int m, double norm)
{

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            c[i * m + j] = 0.0;
            if (fabs(a[i * m + j]) < 1e-50 * norm)
            {
                a[i * m + j] = 0.;
            }
            if (fabs(b[i * m + j]) < 1e-50 * norm)
            {
                b[i * m + j] = 0.;
            }
        }
    }

    int t, q, r;
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    int v3 = m % 3;
    int h3 = m % 3;
    for (r = 0; r < v3; r++) // номер строки
    {
        for (t = 0; t < h3; t++) // номер столбца
        {
            s00 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
            }
            c[r * m + t] = s00;
        }
        for (; t < m; t += 3)
        {
            s00 = 0;
            s01 = 0;
            s02 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
                s01 += a[r * m + q] * b[q * m + t + 1];
                s02 += a[r * m + q] * b[q * m + t + 2];
            }
            c[r * m + t] = s00;
            c[r * m + t + 1] = s01;
            c[r * m + t + 2] = s02;
        }
    }
    for (; r < m; r += 3)
    {
        for (t = 0; t < h3; t++)
        {
            s00 = 0;
            s10 = 0;
            s20 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
                s10 += a[(r + 1) * m + q] * b[q * m + t];
                s20 += a[(r + 2) * m + q] * b[q * m + t];
            }
            c[r * m + t] = s00;
            c[(r + 1) * m + t] = s10;
            c[(r + 2) * m + t] = s20;
        }
        for (; t < m; t += 3)
        {
            s00 = 0;
            s01 = 0;
            s02 = 0;
            s10 = 0;
            s11 = 0;
            s12 = 0;
            s20 = 0;
            s21 = 0;
            s22 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
                s01 += a[r * m + q] * b[q * m + t + 1];
                s02 += a[r * m + q] * b[q * m + t + 2];
                s10 += a[(r + 1) * m + q] * b[q * m + t];
                s11 += a[(r + 1) * m + q] * b[q * m + t + 1];
                s12 += a[(r + 1) * m + q] * b[q * m + t + 2];
                s20 += a[(r + 2) * m + q] * b[q * m + t];
                s21 += a[(r + 2) * m + q] * b[q * m + t + 1];
                s22 += a[(r + 2) * m + q] * b[q * m + t + 2];
            }
            c[r * m + t] = s00;
            c[r * m + t + 1] = s01;
            c[r * m + t + 2] = s02;
            c[(r + 1) * m + t] = s10;
            c[(r + 1) * m + t + 1] = s11;
            c[(r + 1) * m + t + 2] = s12;
            c[(r + 2) * m + t] = s20;
            c[(r + 2) * m + t + 1] = s21;
            c[(r + 2) * m + t + 2] = s22;
        }
    }
}

void matrix_diff(double *a, double *b, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            a[i * col + j] -= b[i * col + j];
        }
    }
}

#include <cmath>

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

int solve_function_task_14(double *matrix, double *inverse_matrix, double *block,
                           double *inv_block, double *help_block, int *perest,
                           int n, int m, double norm_of_matrix)
{
    int k = n / m;                 // количество блоков m на m
    int l = n - k * m;             // размер для прямоугольных блоков
    int bl = (l != 0 ? k + 1 : k); // предел индексирования
    int i_min;                     // индекс блока с минимальной обратной
    int j_min;                     // индекс блока с минимальной обратной

    // инициализируем строку perest

    for (int i = 0; i < k; i++)
    {
        perest[i] = i;
    }

    // printf("\n Matrix \n");
    // printMatrix(matrix, n, n, n);

    for (int p = 0; p < bl; p++) // p - номер шага алгоритма
    {
        double min_norm = __DBL_MAX__; // минимальная норма обратного блока
        int count = 0;                 // счетчик для блоков с обратной матрицей
        i_min = p;
        j_min = p;

        // Находим блок с минимальной нормой обратной и сдвигаем его в верхний левый угол подматрицы
        for (int i = p; i < k; i++) // (i , j) - индекс блока
        {
            for (int j = p; j < k; j++) // размер C_ij = v на h
            {
                get_block(matrix, block, n, m, i, j);
                if (inverse_matrixx(block, block, m, m, norm_of_matrix))
                {
                    // printf("block after inverse: \n");
                    // printMatrix(block, m, m, m);
                    count++;

                    double buf = matrix_norm(block, m, m);
                    if (buf < min_norm)
                    {
                        min_norm = buf;
                        i_min = i;
                        j_min = j;
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
        }

        if (p == k)
        {
            get_block(matrix, block, n, m, p, p);
            if (inverse_matrixx(block, block, l, m, norm_of_matrix))
            {
                // printf("block after inverse: \n");
                // printMatrix(block, m, m, m);
                count++;
                i_min = p;
                j_min = p;
                for (int ii = 0; ii < l; ii++)
                {
                    for (int jj = 0; jj < l; jj++)
                    {
                        inv_block[ii * m + jj] = block[ii * m + jj]; // обратный блок с минимальной нормой
                    }
                }
            }
        }

        if (count == 0)
        {
            printf("ERROR in Solve Function: No block with inverse!\n");
            return 1;
        }

        // printf("\n(%d, %d)\n", i_min, j_min);
        // printMatrix(inv_block, m, m, m);
        if (p != k)
        {
            // меняем местами блоки (i_min, j_min) и (p, p), где p номер шага алгоритма
            change_block_place(matrix, n, m, i_min, j_min, p, p);
            // аналогично надо сделать с приписанной матрицей
            change_block_place(inverse_matrix, n, m, i_min, p, p, p);

            int tmp = perest[p];
            perest[p] = perest[j_min];
            perest[j_min] = tmp;
        }

        // домножаем на обратную в нужной строке и потом вычитаем

        for (int i = p; i < bl; i++)
        {
            for (int j = p; j < bl; j++)
            {
                // int v = (i < k ? m : l); // размер блока (i, j)
                // int h = (j < k ? m : l); // размер блока (i, j)
                if (i == p)
                {
                    get_block(matrix, block, n, m, i, j);
                    if (j == p)
                    {
                        for (int ii = 0; ii < m; ii++)
                        {
                            for (int jj = 0; jj < m; jj++)
                            {
                                block[ii * m + jj] = (ii == jj ? 1 : 0);
                            }
                        }

                        put_block(matrix, block, n, m, i, j);

                        for (int jj = 0; jj < bl; jj++)
                        {
                            // int hh = (jj < k ? m : l); // размер блока (i, jj)
                            get_block(inverse_matrix, block, n, m, i, jj);

                            matrix_mult(inv_block, block, help_block, /*v_p, v_p, v, hh,*/ m, norm_of_matrix);

                            put_block(inverse_matrix, help_block, n, m, i, jj);
                        }

                        for (int jj = p + 1; jj < bl; jj++)
                        {
                            // int hh = (jj < k ? m : l); // размер блока (i, jj)
                            get_block(matrix, block, n, m, i, jj);
                            matrix_mult(inv_block, block, help_block, /*v_p, v_p, v, hh,*/ m, norm_of_matrix);
                            put_block(matrix, help_block, n, m, i, jj);
                        }

                        // printf("matrix after sdapo\n");
                        // printMatrix(matrix, n, n, n);
                        // printf("\nINverse matrix after sdapo\n");
                        // printMatrix(inverse_matrix, n, n, n);

                        continue;
                    }
                }
                else
                {

                    // Вычитаем нашу строку из строк ниже.

                    if (j == p)
                    {
                        get_block(matrix, inv_block, n, m, i, p);
                        for (int jj = 0; jj < bl; jj++)
                        {
                            get_block(inverse_matrix, help_block, n, m, p, jj);
                            matrix_mult(inv_block, help_block, block, /*v, v_p, m, m,*/ m, norm_of_matrix);
                            get_block(inverse_matrix, help_block, n, m, i, jj);
                            matrix_diff(help_block, block, m, m);
                            put_block(inverse_matrix, help_block, n, m, i, jj);
                        }
                    }

                    get_block(matrix, help_block, n, m, p, j);
                    matrix_mult(inv_block, help_block, block, /*v, v_p, m, m,*/ m, norm_of_matrix);
                    get_block(matrix, help_block, n, m, i, j);
                    matrix_diff(help_block, block, m, m);
                    put_block(matrix, help_block, n, m, i, j);
                }
            }
        }

        // printf("matrix after vichet\n");
        // printMatrix(matrix, n, n, n);
        // printf("\nINverse matrix after vichet\n");
        // printMatrix(inverse_matrix, n, n, n);
    }

    // обратный ход
    for (int i = bl - 1; i > 0; i--)
    {
        // int v = (i < k ? m : l); // размер блока (i, j)

        for (int j = i - 1; j >= 0; j--)
        {
            // int b = (j < k ? m : l);
            get_block(matrix, inv_block, n, m, j, i);

            for (int jj = 0; jj < bl; jj++)
            {
                // int h = (jj < k ? m : l); // размер блока (i, jj)

                get_block(inverse_matrix, help_block, n, m, i, jj);
                matrix_mult(inv_block, help_block, block, /*b, v, v, h,*/ m, norm_of_matrix);
                get_block(inverse_matrix, help_block, n, m, j, jj);
                matrix_diff(help_block, block, m, m);
                put_block(inverse_matrix, help_block, n, m, j, jj);
            }
        }
    }

    for (int i = 0; i < k; i++)
    {
        for (int j = i; j < k; j++)
        {
            if (perest[j] == i)
            {
                change_block_place(inverse_matrix, n, m, i, i, j, i);
                int tmp = perest[i];
                perest[i] = perest[j];
                perest[j] = tmp;
            }
        }
    }

    // printf("matrix after obrant\n");
    // printMatrix(matrix, n, n, n);
    // printf("\nINverse matrix after brant\n");
    // printMatrix(inverse_matrix, n, n, n);

    return 0;
}