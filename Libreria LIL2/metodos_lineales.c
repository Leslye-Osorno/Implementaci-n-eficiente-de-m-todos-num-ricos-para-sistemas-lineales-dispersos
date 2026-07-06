#include "metodos_lineales.h"
#include "tools.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

double g_memoria_metodo_kb = 0.0;

void metodoPLU(MatrixLIL* A_lil, double* b, double* x, int* iter, double* tiempo) {
    clock_t inicio = clock();
    int n = A_lil->n;

    // HPC: Estimación de memoria para evitar realloc en el bucle
    int nnz_estimado = A_lil->nnz * 3; // Heurística típica de Fill-in
    if (nnz_estimado < n * 2) nnz_estimado = n * 5;
    
    int cap_L = nnz_estimado, cap_U = nnz_estimado;
    MatrixCSR* L = (MatrixCSR*)malloc(sizeof(MatrixCSR));
    MatrixCSR* U = (MatrixCSR*)malloc(sizeof(MatrixCSR));
    
    L->n = n; L->m = n; U->n = n; U->m = n;
    L->valores = (double*)malloc(cap_L * sizeof(double)); U->valores = (double*)malloc(cap_U * sizeof(double));
    L->columnas = (int*)malloc(cap_L * sizeof(int));      U->columnas = (int*)malloc(cap_U * sizeof(int));
    L->apuntador_filas = (int*)calloc(n + 1, sizeof(int));U->apuntador_filas = (int*)calloc(n + 1, sizeof(int));
    L->inv_diag = (double*)calloc(n, sizeof(double));     U->inv_diag = (double*)calloc(n, sizeof(double));
    
    int nnz_L = 0, nnz_U = 0;
    double* w = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < A_lil->filas[i].tamano; k++) {
            w[A_lil->filas[i].datos[k].col] = A_lil->filas[i].datos[k].val;
        }

        for (int j = 0; j < i; j++) {
            if (fabs(w[j]) > 1e-15) {
                w[j] *= U->inv_diag[j]; 
                for (int ptr = U->apuntador_filas[j]; ptr < U->apuntador_filas[j+1]; ptr++) {
                    int col = U->columnas[ptr];
                    if (col > j) w[col] -= w[j] * U->valores[ptr];
                }
            }
        }

        L->inv_diag[i] = 1.0; 
        for (int j = 0; j < n; j++) {
            if (fabs(w[j]) > 1e-15 || i == j) {
                if (j < i) { 
                    if (nnz_L >= cap_L) { cap_L = (int)(cap_L * 1.5) + 1; L->valores = (double*)realloc(L->valores, cap_L*sizeof(double)); L->columnas = (int*)realloc(L->columnas, cap_L*sizeof(int)); }
                    L->valores[nnz_L] = w[j]; L->columnas[nnz_L] = j; nnz_L++;
                } else { 
                    if (nnz_U >= cap_U) { cap_U = (int)(cap_U * 1.5) + 1; U->valores = (double*)realloc(U->valores, cap_U*sizeof(double)); U->columnas = (int*)realloc(U->columnas, cap_U*sizeof(int)); }
                    if (j > i) {
                        U->valores[nnz_U] = w[j]; U->columnas[nnz_U] = j; nnz_U++;
                    } else if (i == j) {
                        U->inv_diag[i] = 1.0 / w[j];
                    }
                }
            }
            w[j] = 0.0; 
        }
        L->apuntador_filas[i+1] = nnz_L;
        U->apuntador_filas[i+1] = nnz_U;
    }

    g_memoria_metodo_kb = ((nnz_L + nnz_U) * 12.0 + n * 32.0) / 1024.0;

    double* y = (double*)calloc(n, sizeof(double));
    sustitucionHaciaAdelanteCSR(L, b, y);
    sustitucionHaciaAtrasCSR(U, y, x);

    free(w); free(y);
    liberarMatrizCSR(L); liberarMatrizCSR(U);

    *iter = 1; *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}


void metodoCholesky(MatrixLIL* A_lil, double* b, double* x, int* iter, double* tiempo) {
    clock_t inicio = clock();
    int n = A_lil->n;

    // Estimación para Cholesky (Aprovechamos simetría)
    int cap_L = A_lil->nnz * 2; 
    if (cap_L < n * 2) cap_L = n * 5;

    MatrixCSR* L = (MatrixCSR*)malloc(sizeof(MatrixCSR));
    L->n = n; L->m = n;
    L->valores = (double*)malloc(cap_L * sizeof(double));
    L->columnas = (int*)malloc(cap_L * sizeof(int));
    L->apuntador_filas = (int*)calloc((n + 1), sizeof(int));
    L->inv_diag = (double*)calloc(n, sizeof(double));
    int nnz_L = 0;

    double* w = (double*)calloc(n, sizeof(double)); 

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < A_lil->filas[i].tamano; k++) {
            int col = A_lil->filas[i].datos[k].col;
            if (col <= i) w[col] = A_lil->filas[i].datos[k].val;
        }

        for (int j = 0; j <= i; j++) {
            double suma = 0.0;
            int ptr = L->apuntador_filas[j];
            int fin = L->apuntador_filas[j+1];
            for (; ptr < fin; ptr++) {
                int col_k = L->columnas[ptr];
                if (col_k >= j) break;
                suma += w[col_k] * L->valores[ptr];
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

        for (int j = 0; j < i; j++) { 
            if (fabs(w[j]) > 1e-15) {
                if (nnz_L >= cap_L) {
                    cap_L = (int)(cap_L * 1.5) + 1;
                    L->valores = (double*)realloc(L->valores, cap_L * sizeof(double));
                    L->columnas = (int*)realloc(L->columnas, cap_L * sizeof(int));
                }
                L->valores[nnz_L] = w[j];
                L->columnas[nnz_L] = j;
                nnz_L++;
            }
            w[j] = 0.0;
        }
        w[i] = 0.0; 
        L->apuntador_filas[i+1] = nnz_L;
    }

    g_memoria_metodo_kb = (nnz_L * 12.0 + n * 16.0) / 1024.0;

    double* y = (double*)calloc(n, sizeof(double));
    sustitucionHaciaAdelanteCSR(L, b, y);

    // L^T * x = y usando columnas orientadas de L inferior (Mantiene solución exacta sin trasponer memoria)
    for (int i = 0; i < n; i++) x[i] = y[i];
    for (int i = n - 1; i >= 0; i--) {
        x[i] *= L->inv_diag[i];
        for (int ptr = L->apuntador_filas[i]; ptr < L->apuntador_filas[i+1]; ptr++) {
            int col = L->columnas[ptr];
            if (col < i) x[col] -= L->valores[ptr] * x[i];
        }
    }

    free(w); free(y); liberarMatrizCSR(L);
    *iter = 1; *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}


void metodoGradienteConjugado(MatrixLIL* A_lil, double* b, double* x, int* iter, double* tiempo) {
    clock_t inicio = clock();
    int n = A_lil->n;

    MatrixCSR* A = convertirLILaCSR(A_lil);

    double *r = (double*)calloc(n, sizeof(double));
    double *p = (double*)calloc(n, sizeof(double));
    double *Ap = (double*)calloc(n, sizeof(double));
    double *z = (double*)calloc(n, sizeof(double)); 

    for (int i = 0; i < n; i++) { 
        x[i] = 0.0; 
        r[i] = b[i]; 
        z[i] = r[i] * A->inv_diag[i]; 
        p[i] = z[i]; 
    }

    double rz_old = productoPunto(n, r, z);
    int k;

    for (k = 0; k < n; k++) {
        multiplicarMatrizVectorCSR(A, p, Ap); 
        
        double pAp = productoPunto(n, p, Ap);
        double alpha = rz_old / pAp;

        saxpy(n, alpha, p, x, x);
        saxpy(n, -alpha, Ap, r, r);

        double rsnew = productoPunto(n, r, r);
        if (rsnew < 1e-20) break;

        for (int i = 0; i < n; i++) z[i] = r[i] * A->inv_diag[i]; 

        double rz_new = productoPunto(n, r, z);
        double beta = rz_new / rz_old;

        saxpy(n, beta, p, z, p);
        rz_old = rz_new;
    }

    g_memoria_metodo_kb = (A->nnz * 12.0 + n * 48.0) / 1024.0;
    *iter = k; *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;

    free(r); free(p); free(Ap); free(z); liberarMatrizCSR(A);
}
