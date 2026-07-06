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
    double* valores;      
    int* columnas;        
    int* apuntador_filas; 
    double* inv_diag;     
} MatrixCSR;

// --- Gestión de Matriz LIL ---
MatrixLIL* asignarMatrizLIL(int n, int m);
void liberarMatrizLIL(MatrixLIL* A);
void insertarValorLIL(MatrixLIL* A, int f, int c, double valor);
double obtenerValorLIL(MatrixLIL* A, int row, int col);
void multiplicarMatrizVectorLIL(MatrixLIL* A, double* x, double* res);
void sustitucionHaciaAdelante(int n, MatrixLIL* L, double* b, double* y);
void sustitucionHaciaAtras(int n, MatrixLIL* U, double* y, double* x);

// NUEVAS OPTIMIZACIONES LIL
void ordenarLIL(MatrixLIL* A);
void ajustarMemoriaLIL(MatrixLIL* A);

// --- Gestión de Matriz CSR (Optimizada) ---
MatrixCSR* convertirLILaCSR(MatrixLIL* A);
void liberarMatrizCSR(MatrixCSR* csr);
void multiplicarMatrizVectorCSR(MatrixCSR* A, double* x, double* res);
void sustitucionHaciaAdelanteCSR(MatrixCSR* L, double* b, double* y);
void sustitucionHaciaAtrasCSR(MatrixCSR* U, double* y, double* x);

// --- Álgebra Vectorial (SIMD Vectorizado con restrict) ---
void sumarVectores(int n, double* restrict v1, double* restrict v2, double* restrict res);
void restarVectores(int n, double* restrict v1, double* restrict v2, double* restrict res);
double productoPunto(int n, double* restrict v1, double* restrict v2);
double norma2(int n, double* restrict v);
void saxpy(int n, double a, double* restrict x, double* restrict y, double* restrict res);

#endif
