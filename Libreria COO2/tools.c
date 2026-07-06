#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "operaciones_basicas.h"

double calcularResidualFinal(MatrixCOO* A_coo, double* b, double* x) {
    int n = A_coo->n;
    double* Ax = (double*)calloc(n, sizeof(double));
    
    if (Ax == NULL) return -1.0; 

    // Usamos el puente temporal para multiplicar a máxima velocidad
    MatrixCSR* A_csr = convertirCOOaCSR(A_coo);
    multiplicarMatrizVectorCSR(A_csr, x, Ax);
    
    double suma = 0;
    for (int i = 0; i < n; i++) {
        double r = b[i] - Ax[i];
        suma += r * r;
    }
    
    free(Ax);
    liberarMatrizCSR(A_csr);
    return sqrt(suma);
}
