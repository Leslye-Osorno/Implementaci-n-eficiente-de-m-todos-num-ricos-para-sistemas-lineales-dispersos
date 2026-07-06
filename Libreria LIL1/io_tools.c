#include "io_tools.h"
#include <stdio.h>

void registrarCSV(const char* m, int n, int nnz, double mem, double dens, double tL, double tM, double res) {
    FILE* f = fopen("estadisticos_libreriaLIL.csv", "a");
    if (f) {
        fprintf(f, "%s,%d,%d,%.2f,%.6f,%.6f,%.6f,%.2e\n", m, n, nnz, mem, dens, tL, tM, res);
        fclose(f);
    }
}

void guardarTxtDetallado(const char* metodo, const char* matrizNom, int n, double* x) {
    char nom[512];
    sprintf(nom, "sol_%s_%s", metodo, matrizNom);
    FILE* f = fopen(nom, "w");
    if (f) {
        for (int i = 0; i < n; i++) fprintf(f, "%.15f\n", x[i]);
        fclose(f);
    }
}
