#ifndef METODOS_LINEALES_H
#define METODOS_LINEALES_H
#include "operaciones_basicas.h"

extern double g_memoria_metodo_kb;

void metodoPLU(MatrixLIL* A, double* b, double* x, int* iter, double* tiempo);
void metodoCholesky(MatrixLIL* A, double* b, double* x, int* iter, double* tiempo);
void metodoGradienteConjugado(MatrixLIL* A, double* b, double* x, int* iter, double* tiempo);

#endif
