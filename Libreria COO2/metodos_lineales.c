#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "metodos_lineales.h"
#include "operaciones_basicas.h"
#include "errores.h"

void metodoPLU(MatrixCOO* A_coo, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada) {
    if (errorMatrizNoCuadrada(A_coo)) return;
    clock_t inicio = clock();
    int n = A_coo->n;

    //Convertir de COO a CSR 
    MatrixCSR* A = convertirCOOaCSR(A_coo);

    // Estructuras CSR dinámicas para L y U (Crecen solo lo necesario)
    int cap_L = n * 5, cap_U = n * 5; 
    MatrixCSR* L = (MatrixCSR*)malloc(sizeof(MatrixCSR)); L->n = n; L->m = n;
    MatrixCSR* U = (MatrixCSR*)malloc(sizeof(MatrixCSR)); U->n = n; U->m = n;
    
    L->valores = (double*)malloc(cap_L * sizeof(double)); L->columnas = (int*)malloc(cap_L * sizeof(int));
    U->valores = (double*)malloc(cap_U * sizeof(double)); U->columnas = (int*)malloc(cap_U * sizeof(int));
    L->apuntador_filas = (int*)calloc(n + 1, sizeof(int)); U->apuntador_filas = (int*)calloc(n + 1, sizeof(int));
    L->inv_diag = (double*)calloc(n, sizeof(double));      U->inv_diag = (double*)calloc(n, sizeof(double));
    
    int nnz_L = 0, nnz_U = 0;

    // EL ACUMULADOR DENSO (SPA) - La magia del HPC
    double* w = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        //SCATTER (Esparcir la fila de A al acumulador)
        for (int k = A->apuntador_filas[i]; k < A->apuntador_filas[i+1]; k++) {
            w[A->columnas[k]] = A->valores[k];
        }

        //COMPUTE (Operaciones Doolittle en O(1))
        for (int j = 0; j < i; j++) {
            if (fabs(w[j]) > 1e-15) {
                w[j] *= U->inv_diag[j]; // División convertida a multiplicación
                for (int ptr = U->apuntador_filas[j]; ptr < U->apuntador_filas[j+1]; ptr++) {
                    int col = U->columnas[ptr];
                    if (col > j) w[col] -= w[j] * U->valores[ptr]; // Fill-in automático
                }
            }
        }

        // GATHER (Recolectar y limpiar)
        L->inv_diag[i] = 1.0; // Diagonal de L implícita
        for (int j = 0; j < n; j++) {
            if (fabs(w[j]) > 1e-15 || i == j) {
                if (j < i) { // Va para L
                    if (nnz_L >= cap_L) { cap_L *= 2; L->valores = realloc(L->valores, cap_L*sizeof(double)); L->columnas = realloc(L->columnas, cap_L*sizeof(int)); }
                    L->valores[nnz_L] = w[j]; L->columnas[nnz_L] = j; nnz_L++;
                } else { // Va para U
                    if (nnz_U >= cap_U) { cap_U *= 2; U->valores = realloc(U->valores, cap_U*sizeof(double)); U->columnas = realloc(U->columnas, cap_U*sizeof(int)); }
                    if (j > i) {
                        U->valores[nnz_U] = w[j]; U->columnas[nnz_U] = j; nnz_U++;
                    } else if (i == j) {
                        U->inv_diag[i] = 1.0 / w[j]; // Extraer diagonal
                    }
                }
            }
            w[j] = 0.0; // Resetear acumulador
        }
        L->apuntador_filas[i+1] = nnz_L;
        U->apuntador_filas[i+1] = nnz_U;
    }
    L->nnz = nnz_L; U->nnz = nnz_U;
    
    // Guardar memoria usada en el puntero proporcionado
    *memoriaUsada = ((nnz_L + nnz_U) * 12.0 + n * 32.0) / 1024.0;

    // RESOLUCIÓN (Branchless)
    double* y = (double*)calloc(n, sizeof(double));
    sustitucionHaciaAdelanteCSR(L, b, y);
    sustitucionHaciaAtrasCSR(U, y, x);

    free(w); free(y);
    liberarMatrizCSR(A); liberarMatrizCSR(L); liberarMatrizCSR(U);

    *iter = 1; *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}


void metodoCholesky(MatrixCOO* A_coo, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada) {
    if (errorMatrizNoCuadrada(A_coo)) return;
    clock_t inicio = clock();
    int n = A_coo->n;

    MatrixCSR* A = convertirCOOaCSR(A_coo);

    int cap_L = n * 5;
    MatrixCSR* L = (MatrixCSR*)malloc(sizeof(MatrixCSR)); L->n = n; L->m = n;
    L->valores = (double*)malloc(cap_L * sizeof(double));
    L->columnas = (int*)malloc(cap_L * sizeof(int));
    L->apuntador_filas = (int*)calloc((n + 1), sizeof(int));
    L->inv_diag = (double*)calloc(n, sizeof(double));
    int nnz_L = 0;

    double* w = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        // SCATTER
        for (int k = A->apuntador_filas[i]; k < A->apuntador_filas[i+1]; k++) {
            if (A->columnas[k] <= i) w[A->columnas[k]] = A->valores[k];
        }

        // COMPUTE
        for (int j = 0; j <= i; j++) {
            double suma = 0.0;
            for (int ptr = L->apuntador_filas[j]; ptr < L->apuntador_filas[j+1]; ptr++) {
                int c = L->columnas[ptr];
                if (c >= j) break;
                suma += w[c] * L->valores[ptr];
            }

            double val = w[j] - suma;
            if (i == j) {
                double diag_val = sqrt(val);
                w[i] = diag_val;
                L->inv_diag[i] = 1.0 / diag_val;
            } else {
                w[j] = val * L->inv_diag[j];
            }
        }

        // GATHER
        for (int j = 0; j < i; j++) { 
            if (fabs(w[j]) > 1e-15) {
                if (nnz_L >= cap_L) { cap_L *= 2; L->valores = realloc(L->valores, cap_L * sizeof(double)); L->columnas = realloc(L->columnas, cap_L * sizeof(int)); }
                L->valores[nnz_L] = w[j]; L->columnas[nnz_L] = j; nnz_L++;
            }
            w[j] = 0.0;
        }
        w[i] = 0.0; 
        L->apuntador_filas[i+1] = nnz_L;
    }
    L->nnz = nnz_L;
    
    // Guardar memoria usada en el puntero proporcionado
    *memoriaUsada = (nnz_L * 12.0 + n * 16.0) / 1024.0;

    // RESOLUCIÓN (Sustitución + Transpuesta Implícita con Saxpy)
    double* y = (double*)calloc(n, sizeof(double));
    sustitucionHaciaAdelanteCSR(L, b, y);

    for (int i = 0; i < n; i++) x[i] = y[i];
    for (int i = n - 1; i >= 0; i--) {
        x[i] *= L->inv_diag[i];
        for (int ptr = L->apuntador_filas[i]; ptr < L->apuntador_filas[i+1]; ptr++) {
            int col = L->columnas[ptr];
            if (col < i) x[col] -= L->valores[ptr] * x[i];
        }
    }

    free(w); free(y); liberarMatrizCSR(A); liberarMatrizCSR(L);
    *iter = 1; *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}


void gradienteConjugado(MatrixCOO* A_coo, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada) {
    clock_t inicio = clock();
    int n = A_coo->n;

    // Puente a Caché y Extracción del Precondicionador
    MatrixCSR* A = convertirCOOaCSR(A_coo);

    double *r = (double*)calloc(n, sizeof(double));
    double *p = (double*)calloc(n, sizeof(double));
    double *Ap = (double*)calloc(n, sizeof(double));
    double *z = (double*)calloc(n, sizeof(double)); 

    for (int i = 0; i < n; i++) { 
        x[i] = 0.0; 
        r[i] = b[i]; 
        z[i] = r[i] * A->inv_diag[i]; // Precondicionador Inicial
        p[i] = z[i]; 
    }

    double rz_old = productoPuntoDisperso(n, r, z);
    int k_iter;

    for (k_iter = 0; k_iter < n; k_iter++) {
        // Multiplicación en Caché L1 (Rapidísima)
        multiplicarMatrizVectorCSR(A, p, Ap);
        
        double pAp = productoPuntoDisperso(n, p, Ap);
        double alpha = rz_old / pAp;

        saxpy(n, alpha, p, x, x);
        saxpy(n, -alpha, Ap, r, r);

        double rsnew = productoPuntoDisperso(n, r, r);
        if (rsnew < 1e-20) break; // Adiós sqrt()

        // Aplicar Precondicionador de Jacobi (Aplasta iteraciones)
        for (int i = 0; i < n; i++) z[i] = r[i] * A->inv_diag[i];

        double rz_new = productoPuntoDisperso(n, r, z);
        double beta = rz_new / rz_old;

        saxpy(n, beta, p, z, p);
        rz_old = rz_new;
    }

    // Guardar memoria usada en el puntero proporcionado
    *memoriaUsada = (A->nnz * 12.0 + n * 48.0) / 1024.0;
    
    *iter = k_iter;
    *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;

    free(r); free(p); free(Ap); free(z); liberarMatrizCSR(A);
}
