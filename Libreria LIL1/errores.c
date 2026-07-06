#include "errores.h"

int errorMatrizNoCuadrada(MatrixLIL* A) { return (A->n != A->m); }

int verificarMatrizSimetrica(MatrixLIL* A) {
    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->filas[i].tamano; j++) {
            int col = A->filas[i].datos[j].col;
            double val = A->filas[i].datos[j].val;
            if (fabs(val - obtenerValorLIL(A, col, i)) > 1e-12) return 1;
        }
    }
    return 0;
}
