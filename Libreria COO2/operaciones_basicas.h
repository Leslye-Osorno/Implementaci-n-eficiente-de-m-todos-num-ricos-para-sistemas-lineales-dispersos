#ifndef OPERACIONES_BASICAS_H
#define OPERACIONES_BASICAS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Formato de Ingesta Rápida (Lectura) ---
typedef struct {
    int n, m, nnz;
    int* row;
    int* col;
    double* val;
    int last_idx;
} MatrixCOO;

// --- NUEVO: Formato HPC Matemático ---
typedef struct {
    int n, m, nnz;
    double* valores;      // Memoria 100% contigua
    int* columnas;        // Índices contiguos
    int* apuntador_filas; // Punteros de inicio y fin de fila
    double* inv_diag;     // Precondicionador Jacobi y eliminador de "Ifs"
} MatrixCSR;

// Gestión COO original
MatrixCOO* asignarMatrizCOO(int n, int m, int nnz);
void liberarMatrizCOO(MatrixCOO* A);
void inicializarMatrizDispersa(MatrixCOO* matrix);
double obtenerValorCOO(MatrixCOO* A, int row, int col);

// --- NUEVO: El Puente a HPC ---
MatrixCSR* convertirCOOaCSR(MatrixCOO* A);
void liberarMatrizCSR(MatrixCSR* csr);

// --- Álgebra Dispersa Optimizada (Reemplaza a las listas enlazadas) ---
void multiplicarMatrizVectorCSR(MatrixCSR* A, double* x, double* res);
void sustitucionHaciaAdelanteCSR(MatrixCSR* L, double* b, double* y);
void sustitucionHaciaAtrasCSR(MatrixCSR* U, double* y, double* x);

// Álgebra Vectorial (Tus funciones originales intáctas)
double productoPuntoDisperso(int n, double* v1, double* v2);
void multiplicarVectorEscalar(int n, double escalar, double* v, double* resultado);
// (Si necesitas saxpy o sumarVectores, puedes mantenerlas aquí)
void saxpy(int n, double a, double* x, double* y, double* res);

#endif
