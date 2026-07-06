#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "leer_matriz.h"
#include "operaciones_basicas.h"
#include "errores.h"

MatrixCOO* leerMatrizArchivo(const char* archivo, double* tiempoLectura) {
    clock_t inicio = clock();
    FILE* f = fopen(archivo, "r");
    if (f == NULL) {
        printf("Error: No se pudo abrir el archivo %s\n", archivo);
        exit(1);
    }

    int nnz, n, m;
    // Ahora leemos el formato nnz, n, m
    if (fscanf(f, "%d %d %d", &nnz, &n, &m) != 3) {
        printf("Error leyendo dimensiones y nnz (Formato esperado: nnz n m)\n");
        fclose(f);
        exit(1);
    }

    MatrixCOO* A = asignarMatrizCOO(n, m, nnz);

    int k = 0;
    double valor;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (fscanf(f, "%lf", &valor) != 1) {
                printf("Error leyendo elemento [%d,%d]\n", i, j);
                fclose(f);
                exit(1);
            }
            // Solo guardamos valores distintos de cero (formato disperso)
            if (valor != 0.0) {
                A->row[k] = i;
                A->col[k] = j;
                A->val[k] = valor;
                k++;
            }
        }
    }

    fclose(f);
    *tiempoLectura = (double)(clock() - inicio) / CLOCKS_PER_SEC;
    return A;
}
