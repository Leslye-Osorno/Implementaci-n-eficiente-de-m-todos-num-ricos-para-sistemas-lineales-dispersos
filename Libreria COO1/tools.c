#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "operaciones_basicas.h"

double calcularResidualFinal(MatrixCOO* A, double* b, double* x) {
    int n = A->n;
    double* Ax = (double*)calloc(n, sizeof(double));
    
    // Verificamos que la memoria se haya asignado correctamente
    if (Ax == NULL) {
        return -1.0; 
    }

    multiplicarMatrizVectorDisperso(A, x, Ax);
    
    double suma = 0;
    for (int i = 0; i < n; i++) {
        double r = b[i] - Ax[i];
        suma += r * r;
    }
    
    free(Ax);
    return sqrt(suma);
}
