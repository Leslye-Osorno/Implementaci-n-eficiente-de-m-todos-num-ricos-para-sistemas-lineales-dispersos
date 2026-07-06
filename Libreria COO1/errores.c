#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "errores.h"
#include "operaciones_basicas.h"

int errorMatrizNoCuadrada(MatrixCOO* A) {
    if (A->n != A->m) {
        printf("ERROR: La matriz no es cuadrada (%dx%d)\n", A->n, A->m);
        return 1;
    }
    return 0;
}

int errorDiagonalCeros(MatrixCOO* A) {
    for (int i = 0; i < A->n; i++) {
        double d = obtenerValorCOO(A, i, i);
        if (fabs(d) < 1e-15) {
            printf("ERROR: Elemento diagonal cero o muy pequeño en fila %d\n", i);
            return 1;
        }
    }
    return 0;
}

int verificarMatrizSimetrica(MatrixCOO* A) {
    for (int k = 0; k < A->nnz; k++) {
        int i = A->row[k];
        int j = A->col[k];
        double val_ij = A->val[k];
        double val_ji = obtenerValorCOO(A, j, i);

        if (fabs(val_ij - val_ji) > 1e-9) {
            printf("ERROR: La matriz no es simétrica en (%d,%d)\n", i, j);
            return 1;
        }
    }
    return 0;
}

int verificarDefinidaPositiva(MatrixCOO* A) {
    for (int i = 0; i < A->n; i++) {
        double d = obtenerValorCOO(A, i, i);
        if (d <= 0) {
            printf("ADVERTENCIA: La matriz puede no ser positiva definida (diagonal[%d] = %.4f)\n", i, d);
            return 1;
        }
    }
    return 0;
}

int errorDimensionVector(int n, int m) {
    if (n != m) {
        printf("ERROR: Dimensiones incompatibles entre matriz y vector\n");
        return 1;
    }
    return 0;
}
