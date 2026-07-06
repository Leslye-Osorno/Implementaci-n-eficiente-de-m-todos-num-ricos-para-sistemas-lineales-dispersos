#include "leer_matriz.h"
#include <time.h>

MatrixLIL* leerMatrizArchivo(const char* archivo, double* tiempoLectura) {
    clock_t inicio = clock();
    FILE* f = fopen(archivo, "r");
    if (!f) return NULL;
    
    int nnz_ref, n, m;
    if (fscanf(f, "%d %d %d", &nnz_ref, &n, &m) != 3) { fclose(f); return NULL; }
    
    MatrixLIL* A = asignarMatrizLIL(n, m);
    double val;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (fscanf(f, "%lf", &val) == 1 && fabs(val) > 1e-18) {
                insertarValorLIL(A, i, j, val);
            }
        }
    }
    fclose(f);
    
    // Pasos HPC de consolidación de la estructura LIL
    ordenarLIL(A);        // Ordena y agrupa eliminando repeticiones
    ajustarMemoriaLIL(A); // Devuelve al sistema toda la memoria RAM que sobra
    
    *tiempoLectura = (double)(clock() - inicio) / CLOCKS_PER_SEC;
    return A;
}
