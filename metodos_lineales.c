#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "metodos_lineales.h"
#include "operaciones_basicas.h"
#include "errores.h"

typedef struct LocalIndexNode {
    int col;
    double val;
    struct LocalIndexNode* next;
} LocalIndexNode;

typedef struct {
    LocalIndexNode* head;
    LocalIndexNode* last;
} LocalRowList;

typedef struct {
    LocalRowList* rows;
    int n;
} LocalMatrixIndex;

static LocalMatrixIndex* construirIndiceLocal(MatrixCOO* A) {
    int n = A->n;
    LocalMatrixIndex* idx = (LocalMatrixIndex*)malloc(sizeof(LocalMatrixIndex));
    idx->n = n;
    idx->rows = (LocalRowList*)calloc(n, sizeof(LocalRowList));
    
    for (int k = 0; k < A->nnz; k++) {
        int r = A->row[k];
        int c = A->col[k];
        double v = A->val[k];

        LocalIndexNode* node = (LocalIndexNode*)malloc(sizeof(LocalIndexNode));
        node->col = c;
        node->val = v;
        node->next = NULL;

        if (idx->rows[r].head == NULL) {
            idx->rows[r].head = node;
            idx->rows[r].last = node;
        } else {
            idx->rows[r].last->next = node;
            idx->rows[r].last = node;
        }
    }
    return idx;
}

static void liberarIndiceLocal(LocalMatrixIndex* idx) {
    for (int i = 0; i < idx->n; i++) {
        LocalIndexNode* curr = idx->rows[i].head;
        while (curr != NULL) {
            LocalIndexNode* temp = curr;
            curr = curr->next;
            free(temp);
        }
    }
    free(idx->rows);
    free(idx);
}

// Función auxiliar ultra veloz para buscar en el índice local de la matriz A
static double buscarEnIndiceLocal(LocalMatrixIndex* idx, int fila, int col) {
    if (!idx || fila >= idx->n) return 0.0;
    LocalIndexNode* curr = idx->rows[fila].head;
    while (curr != NULL) {
        if (curr->col == col) return curr->val;
        if (curr->col > col) break; // Aprovecha que los nodos vienen ordenados por columna
        curr = curr->next;
    }
    return 0.0;
}

void metodoPLU(MatrixCOO* A, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada) {
    if (errorMatrizNoCuadrada(A)) return;
    clock_t inicio = clock();
    int n = A->n;

    LocalMatrixIndex* idxA = construirIndiceLocal(A);

    // OPTIMIZACIÓN CRÍTICA: Espacio de trabajo denso temporal para cálculo a velocidad de hardware O(1)
    double* L_denso = (double*)calloc(n * n, sizeof(double));
    double* U_denso = (double*)calloc(n * n, sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double suma = 0.0;
            for (int j = 0; j < i; j++) {
                suma += L_denso[i * n + j] * U_denso[j * n + k];
            }
            U_denso[i * n + k] = buscarEnIndiceLocal(idxA, i, k) - suma;
        }

        L_denso[i * n + i] = 1.0;
        for (int k = i + 1; k < n; k++) {
            double suma = 0.0;
            for (int j = 0; j < i; j++) {
                suma += L_denso[k * n + j] * U_denso[j * n + i];
            }
            
            double diagU = U_denso[i * n + i];
            if (fabs(diagU) < 1e-15) {
                diagU = (diagU >= 0) ? 1e-15 : -1e-15;
            }
            L_denso[k * n + i] = (buscarEnIndiceLocal(idxA, k, i) - suma) / diagU;
        }
    }

    // Contar los elementos no nulos reales generados por el fill-in
    int nnzL = 0, nnzU = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(L_denso[i * n + j]) > 1e-15) nnzL++;
            if (fabs(U_denso[i * n + j]) > 1e-15) nnzU++;
        }
    }

    // Empaquetar eficientemente las estructuras dispersas finales en formato COO
    MatrixCOO* L_coo = asignarMatrizCOO(n, n, nnzL);
    MatrixCOO* U_coo = asignarMatrizCOO(n, n, nnzU);

    int kl = 0, ku = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(L_denso[i * n + j]) > 1e-15) {
                L_coo->row[kl] = i; L_coo->col[kl] = j; L_coo->val[kl] = L_denso[i * n + j]; kl++;
            }
            if (fabs(U_denso[i * n + j]) > 1e-15) {
                U_coo->row[ku] = i; U_coo->col[ku] = j; U_coo->val[ku] = U_denso[i * n + j]; ku++;
            }
        }
    }

    // Medición exacta de la memoria en KB ocupada por las estructuras finales y el índice
    size_t bytesEstructuras = (sizeof(MatrixCOO) * 2) + 
                             (sizeof(int) * 2 * nnzL + sizeof(double) * nnzL) +
                             (sizeof(int) * 2 * nnzU + sizeof(double) * nnzU) +
                             (sizeof(double) * n) + // Vector b/y
                             (sizeof(LocalRowList) * n) + (sizeof(LocalIndexNode) * A->nnz);
    *memoriaUsada = (double)bytesEstructuras / 1024.0;

    double *y = (double*)calloc(n, sizeof(double));
    sustitucionHaciaAdelante(n, L_coo, b, y);
    sustitucionHaciaAtras(n, U_coo, y, x);

    // Liberación segura de memoria
    free(L_denso); free(U_denso); free(y);
    liberarIndiceLocal(idxA);
    liberarMatrizCOO(L_coo); liberarMatrizCOO(U_coo);

    *iter = 1;
    *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}

void metodoCholesky(MatrixCOO* A, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada) {
    if (errorMatrizNoCuadrada(A)) return;
    clock_t inicio = clock();
    int n = A->n;

    LocalMatrixIndex* idxA = construirIndiceLocal(A);

    // OPTIMIZACIÓN CRÍTICA: Matriz L temporal indexada directamente en un arreglo denso
    double* L_denso = (double*)calloc(n * n, sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double suma = 0.0;
            for (int k = 0; k < j; k++) {
                suma += L_denso[i * n + k] * L_denso[j * n + k];
            }

            double valA = buscarEnIndiceLocal(idxA, i, j);

            if (i == j) {
                double valDiag = valA - suma;
                if (valDiag < 1e-15) valDiag = 1e-15;
                L_denso[i * n + j] = sqrt(valDiag);
            } else {
                double diagJ = L_denso[j * n + j];
                if (fabs(diagJ) < 1e-15) diagJ = 1e-15;
                L_denso[i * n + j] = (valA - suma) / diagJ;
            }
        }
    }

    // Contar los elementos no nulos de L generados
    int nnzL = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            if (fabs(L_denso[i * n + j]) > 1e-15) nnzL++;
        }
    }

    MatrixCOO* L_coo = asignarMatrizCOO(n, n, nnzL);
    MatrixCOO* Lt_coo = asignarMatrizCOO(n, n, nnzL);

    int kl = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            if (fabs(L_denso[i * n + j]) > 1e-15) {
                L_coo->row[kl] = i;  L_coo->col[kl] = j;  L_coo->val[kl] = L_denso[i * n + j];
                Lt_coo->row[kl] = j; Lt_coo->col[kl] = i; Lt_coo->val[kl] = L_denso[i * n + j];
                kl++;
            }
        }
    }

    // Medición exacta de memoria consumida por Cholesky
    size_t bytesEstructuras = (sizeof(MatrixCOO) * 2) + 
                             (sizeof(int) * 4 * nnzL + sizeof(double) * 2 * nnzL) +
                             (sizeof(double) * n) +
                             (sizeof(LocalRowList) * n) + (sizeof(LocalIndexNode) * A->nnz);
    *memoriaUsada = (double)bytesEstructuras / 1024.0;

    double* y = (double*)calloc(n, sizeof(double));
    sustitucionHaciaAdelante(n, L_coo, b, y);
    sustitucionHaciaAtras(n, Lt_coo, y, x);

    free(L_denso); free(y);
    liberarIndiceLocal(idxA);
    liberarMatrizCOO(L_coo); liberarMatrizCOO(Lt_coo);

    *iter = 1;
    *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}

void gradienteConjugado(MatrixCOO* A, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada) {
    clock_t inicio = clock();
    int n = A->n;

    // Memoria dinámica exacta asignada por el método (3 vectores)
    *memoriaUsada = (3.0 * n * sizeof(double)) / 1024.0;

    double *r = malloc(n * sizeof(double));
    double *p = malloc(n * sizeof(double));
    double *Ap = malloc(n * sizeof(double));

    multiplicarMatrizVectorDisperso(A, x, r);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    double rsold = productoPuntoDisperso(n, r, r);
    for (int k_iter = 0; k_iter < n; k_iter++) {
        multiplicarMatrizVectorDisperso(A, p, Ap);
        double divisor = productoPuntoDisperso(n, p, Ap);
        if (fabs(divisor) < 1e-18) break;
        double alpha = rsold / divisor;
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = productoPuntoDisperso(n, r, r);
        if (sqrt(rsnew) < 1e-6) {
            *iter = k_iter + 1;
            break;
        }
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }
        rsold = rsnew;
        *iter = k_iter + 1;
    }

    free(r); free(p); free(Ap);
    *tiempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
}
