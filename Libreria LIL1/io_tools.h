#ifndef IO_TOOLS_H
#define IO_TOOLS_H

void registrarCSV(const char* m, int n, int nnz, double mem, double dens, double tL, double tM, double res);
void guardarTxtDetallado(const char* metodo, const char* matrizNom, int n, double* x);

#endif
