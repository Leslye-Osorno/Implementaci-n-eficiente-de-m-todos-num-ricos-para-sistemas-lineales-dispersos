#include "operaciones_basicas.h"
#include "errores.h"
#include <string.h>


MatrixCOO* asignarMatrizCOO(int n, int m, int nnz) {
    MatrixCOO* A = (MatrixCOO*)malloc(sizeof(MatrixCOO));
    A->n = n; A->m = m; A->nnz = nnz;
    A->row = (int*)malloc(nnz * sizeof(int));
    A->col = (int*)malloc(nnz * sizeof(int));
    A->val = (double*)malloc(nnz * sizeof(double));
    A->last_idx = 0;
    return A;
}

void liberarMatrizCOO(MatrixCOO* A) {
    if (!A) return;
    free(A->row); free(A->col); free(A->val); free(A);
}

double obtenerValorCOO(MatrixCOO* A, int row, int col) {
    for (int i = 0; i < A->last_idx; i++) {
        if (A->row[i] == row && A->col[i] == col) return A->val[i];
    }
    return 0.0;
}


// MÓDULO CSR DE ALTO RENDIMIENTO (HPC)

MatrixCSR* convertirCOOaCSR(MatrixCOO* A) {
    if (!A) return NULL;
    MatrixCSR* csr = (MatrixCSR*)malloc(sizeof(MatrixCSR));
    int n = A->n;
    csr->n = n; csr->m = A->m; csr->nnz = A->last_idx;

    csr->valores = (double*)malloc(csr->nnz * sizeof(double));
    csr->columnas = (int*)malloc(csr->nnz * sizeof(int));
    csr->apuntador_filas = (int*)calloc((n + 1), sizeof(int));
    csr->inv_diag = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) csr->inv_diag[i] = 1.0;

    // Contar elementos por fila
    for (int i = 0; i < csr->nnz; i++) {
        csr->apuntador_filas[A->row[i] + 1]++;
    }

    // Suma acumulativa para obtener punteros de inicio
    for (int i = 0; i < n; i++) {
        csr->apuntador_filas[i + 1] += csr->apuntador_filas[i];
    }

    // Copia temporal de apuntadores para insertar elementos ordenadamente
    int* offset = (int*)malloc(n * sizeof(int));
    memcpy(offset, csr->apuntador_filas, n * sizeof(int));

    // Acomodar valores y columnas en arreglos contiguos
    for (int i = 0; i < csr->nnz; i++) {
        int r = A->row[i];
        int c = A->col[i];
        double v = A->val[i];

        int dest = offset[r]++;
        csr->columnas[dest] = c;
        csr->valores[dest] = v;

        // Extracción de diagonal para Precondicionador / Sustituciones
        if (r == c && fabs(v) > 1e-15) {
            csr->inv_diag[r] = 1.0 / v;
        }
    }

    free(offset);
    return csr;
}

void liberarMatrizCSR(MatrixCSR* csr) {
    if (!csr) return;
    free(csr->valores); free(csr->columnas);
    free(csr->apuntador_filas); free(csr->inv_diag); free(csr);
}

// Multiplicación SpMV Optimizada (Caché L1 Friendly)
void multiplicarMatrizVectorCSR(MatrixCSR* A, double* x, double* res) {
    for (int i = 0; i < A->n; i++) {
        double suma = 0.0;
        int inicio = A->apuntador_filas[i];
        int fin = A->apuntador_filas[i+1];
        for (int j = inicio; j < fin; j++) {
            suma += A->valores[j] * x[A->columnas[j]];
        }
        res[i] = suma;
    }
}

// Sustituciones ultrarrápidas Branchless (Cero 'ifs')
void sustitucionHaciaAdelanteCSR(MatrixCSR* L, double* b, double* y) {
    for (int i = 0; i < L->n; i++) {
        double suma = 0.0;
        for (int j = L->apuntador_filas[i]; j < L->apuntador_filas[i+1]; j++) {
            suma += L->valores[j] * y[L->columnas[j]];
        }
        y[i] = (b[i] - suma) * L->inv_diag[i]; // Aritmética pura
    }
}

void sustitucionHaciaAtrasCSR(MatrixCSR* U, double* y, double* x) {
    for (int i = U->n - 1; i >= 0; i--) {
        double suma = 0.0;
        for (int j = U->apuntador_filas[i]; j < U->apuntador_filas[i+1]; j++) {
            suma += U->valores[j] * x[U->columnas[j]];
        }
        x[i] = (y[i] - suma) * U->inv_diag[i];
    }
}

// --- Álgebra Vectorial Básica ---
double productoPuntoDisperso(int n, double* v1, double* v2) {
    double r = 0.0;
    for (int i = 0; i < n; i++) r += v1[i] * v2[i];
    return r;
}

void saxpy(int n, double a, double* x, double* y, double* res) {
    for (int i = 0; i < n; i++) res[i] = a * x[i] + y[i];
}
