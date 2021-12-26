#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>
#pragma warning(disable : 4996)

void Analyze(float** matrix, int size, int row, int col, float** newMatrix) {
   int newRow = 0;
   int newCol = 0;
   for (int line = 0; line < size - 1; line++) {
      if (line == row) {
         newRow = 1;
      }
      newCol = 0;
      for (int column = 0; column < size - 1; column++) {

         if (column == col) {
            newCol = 1;
         }
         newMatrix[line][column] = matrix[line + newRow][column + newCol];
      }
   }
}

float Determinant(float** matrix, int size) {
   float det = 0;
   float step = 1;

   if (size == 1)
      return matrix[0][0]; //сводим до единичной матрицы
   else if (size == 2)
      return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]; //сводим до матрицы 2*2
   else {
      float** newMatrix = (float**)malloc(sizeof(int*) * size - 1);;
      for (int line = 0; line < size - 1; line++) {
         newMatrix[line] = (float*)malloc(sizeof(int) * (size - 1));;
      }
      if (!newMatrix) {
         printf("Malloc error!\n");
         return 0;
      }
      for (int column = 0; column < size; column++) {
         Analyze(matrix, size, 0, column, newMatrix);
         det = det + (step * matrix[0][column] * Determinant(newMatrix, size - 1));
         step = -step;
      }
      for (int line = 0; line < size - 1; line++) {
         free(newMatrix[line]);
      }
      free(newMatrix);
   }
   return det;
}

void Transp(float** A, int size) {
   double s;
   for (int line = 0; line < size; line++)
      for (int column = line + 1; column < size; column++) {
         s = A[line][column];
         A[line][column] = A[column][line];
         A[column][line] = s;
      }
}

void Inversion(float** matrix, int size, float det, float** newmat) {
   for (int line = 0; line < size; line++) {
      for (int column = 0; column < size; column++) {
         float** dopmat = (float**)malloc(sizeof(int*) * (size - 1));
         for (int line = 0; line < size - 1; line++)
            dopmat[line] = (float*)malloc(sizeof(int) * (size - 1));
         if (!dopmat) {
            printf("Malloc error!\n");
            return;
         }
         int ext = 0;
         Analyze(matrix, size, line, column, dopmat);
         ext = Determinant(dopmat, size - 1);
         if ((line + column) % 2 == 0)
            newmat[line][column] = ext / det;
         else
            newmat[line][column] = -ext / det;
         free(dopmat);
      }
   }
   Transp(newmat, size);
   printf("Inverse matrix:\n");
   for (int line = 0; line < size; line++) {
      for (int column = 0; column < size; column++) {
         printf("%.2f ", newmat[line][column]);
      }
      printf("\n");
   }
}

void SolveSlY(int sizeSLY, float** SLY, float* SolveSLY, float* Solve) {
   float detMat = Determinant(SLY, sizeSLY);
   for (int line = 0; line < sizeSLY; ++line) {
      float** reserv = (float**)malloc(sizeof(int*) * sizeSLY);
      for (int o = 0; o < sizeSLY; o++)
         reserv[o] = (float*)malloc(sizeof(int) * sizeSLY);
      if (!reserv) {
         printf("Malloc error!\n");
         return;
      }
      for (int k = 0; k < sizeSLY; ++k)
         for (int t = 0; t < sizeSLY; ++t)
            reserv[k][t] = SLY[k][t];
      for (int column = 0; column < sizeSLY; ++column)
         reserv[column][line] = SolveSLY[column];
      Solve[line] = (Determinant(reserv, sizeSLY)) / detMat;
      printf("x[%d] = %0.2f\n", line, Solve[line]);
      free(reserv);
   }
}

int main() {
   setlocale(LC_ALL, "Russian");
   printf("Введите размер матрицы: ");
   int metre;
   scanf_s("%d", &metre); //размер
   float** arr = (float**)malloc(sizeof(int*) * metre);
   for (int line = 0; line < metre; line++)
      arr[line] = (float*)malloc(sizeof(int) * metre);
   if (!arr)
   {
      printf("Malloc error!\n");
      return 0;
   }
   printf("\n");
   printf("Введите матрицу размером %dx%d: \n", metre, metre);
   for (int line = 0; line < metre; line++) {
      printf("");
      for (int column = 0; column < metre; column++) {
         scanf_s("%f", &arr[line][column]);
      }
   }
   printf("\n");

   printf("Введите свободные члены уравнений\n");
   float* Solution = (float*)malloc(sizeof(int*) * metre);
   if (!Solution)
   {
      printf("Malloc error!\n");
      return 0;
   }
   for (int line = 0; line < metre; line++)
      scanf_s("%f", &Solution[line]);
   printf("\n");
   printf("Матрица:\n");

   for (int line = 0; line < metre; line++) {
      for (int column = 0; column < metre; column++) {
         printf("%.0f ", arr[line][column]);
      }
      printf("\n");
   }
   printf("\n");
   printf("Свободные члены:\n");
   for (int line = 0; line < metre; line++)
   {
      printf("%.2f ", Solution[line]);
   }
   printf("\n");
   float Det = Determinant(arr, metre);
   printf("Определитель матрицы: %.0f \n", Det);
   printf("\n");
   float** newMatrix = (float**)malloc(sizeof(int*) * metre);
   if (!newMatrix) {
      printf("Malloc error!\n");
      return 0;
   }
   for (int line = 0; line < metre; line++) {
      newMatrix[line] = (float*)malloc(sizeof(int) * metre);
   }
   Inversion(arr, metre, Det, newMatrix);
   printf("\n");
   float* Solve = (float*)malloc(sizeof(int*) * metre);
   if (!Solve) {
      printf("Malloc error!\n");
      return 0;
   }
   SolveSlY(metre, arr, Solution, Solve);
   return 0;
}