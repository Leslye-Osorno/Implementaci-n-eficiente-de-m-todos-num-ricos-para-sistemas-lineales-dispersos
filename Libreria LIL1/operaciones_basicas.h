#ifndef OPERACIONES_BASICAS_H
#define OPERACIONES_BASICAS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- ESTRUCTURAS ORIGINALES (Para lectura del .txt) ---
typedef struct {
    int col;
    double val;
} ElementoLIL;

typedef struct {
    ElementoLIL* datos;
    int tamano;
    int capacidad;
} FilaLIL;

typedef struct {
    int n, m, nnz;
    FilaLIL* filas;
} MatrixLIL;

// --- ESTRUCTURAS HPC OPTIMIZADAS (Para matemáticas) ---
typedef struct {
    int n, m, nnz;
    double* valores;      // Memoria 100% contigua
    int* columnas;        // Índices contiguos
    int* apuntador_filas; // Punteros de inicio de fila
    double* inv_diag;     // Precondicionador y Branchless solver
} MatrixCSR;

// --- Gestión de Matriz LIL ---
MatrixLIL* asignarMatrizLIL(int n, int m);
void liberarMatrizLIL(MatrixLIL* A);
void insertarValorLIL(MatrixLIL* A, int f, int c, double valor);
double obtenerValorLIL(MatrixLIL* A, int row, int col);
void multiplicarMatrizVectorLIL(MatrixLIL* A, double* x, double* res);
void sustitucionHaciaAdelante(int n, MatrixLIL* L, double* b, double* y);
void sustitucionHaciaAtras(int n, MatrixLIL* U, double* y, double* x);

// --- Gestión de Matriz CSR (Optimizada) ---
MatrixCSR* convertirLILaCSR(MatrixLIL* A);
void liberarMatrizCSR(MatrixCSR* csr);
void multiplicarMatrizVectorCSR(MatrixCSR* A, double* x, double* res);
void sustitucionHaciaAdelanteCSR(MatrixCSR* L, double* b, double* y);
void sustitucionHaciaAtrasCSR(MatrixCSR* U, double* y, double* x);

// --- Álgebra Vectorial ---
double productoPunto(int n, double* v1, double* v2);
double norma2(int n, double* v);
void saxpy(int n, double a, double* x, double* y, double* res);
void sumarVectores(int n, double* v1, double* v2, double* res);
void restarVectores(int n, double* v1, double* v2, double* res);

#endif
