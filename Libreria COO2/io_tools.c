#include <stdio.h>
#include "io_tools.h"

void guardarTxtDetallado(const char* metodo, const char* matrizNom, int n, double* x) {
    char nombreArchivo[512];
    sprintf(nombreArchivo, "resultados_%s_%s", metodo, matrizNom);
    FILE* f = fopen(nombreArchivo, "w");
    if (!f) return;
    fprintf(f, "Metodo: %s | Matriz: %s\n", metodo, matrizNom);
    for (int i = 0; i < n; i++) {
        fprintf(f, "x[%d] = %.15f\n", i, x[i]);
    }
    fclose(f);
}

void guardarCsvResultados(const char* metodo, const char* matrizNom, int n, double* x) {
    char nombreArchivo[512];
    sprintf(nombreArchivo, "resultados_%s_%s.csv", metodo, matrizNom);
    FILE* f = fopen(nombreArchivo, "w");
    if (!f) return;
    fprintf(f, "Indice,Resultado\n");
    for (int i = 0; i < n; i++) {
        fprintf(f, "%d,%.15f\n", i, x[i]);
    }
    fclose(f);
}
