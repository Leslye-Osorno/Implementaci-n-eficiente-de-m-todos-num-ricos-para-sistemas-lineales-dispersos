#include "tools.h"

double obtenerMemoriaMatrizKB(MatrixLIL* A) {
    if (!A) return 0.0;
    size_t bytes = sizeof(MatrixLIL) + (A->n * sizeof(FilaLIL));
    for (int i = 0; i < A->n; i++) bytes += A->filas[i].capacidad * sizeof(ElementoLIL);
    return (double)bytes / 1024.0;
}

double calcularResidualFinal(MatrixLIL* A, double* b, double* x) {
    int n = A->n;
    double* Ax = (double*)calloc(n, sizeof(double));
    multiplicarMatrizVectorLIL(A, x, Ax);
    double suma = 0;
    for (int i = 0; i < n; i++) {
        double r = b[i] - Ax[i];
        suma += r * r;
    }
    free(Ax);
    return sqrt(suma);
}

double calcularResidualFinalCSR(MatrixCSR* A, double* b, double* x) {
    int n = A->n;
    double* Ax = (double*)calloc(n, sizeof(double));
    multiplicarMatrizVectorCSR(A, x, Ax); // Usar la ultra-rápida
    double suma = 0;
    for (int i = 0; i < n; i++) {
        double r = b[i] - Ax[i];
        suma += r * r;
    }
    free(Ax);
    return sqrt(suma);
}
