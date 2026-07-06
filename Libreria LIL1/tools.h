#ifndef TOOLS_H
#define TOOLS_H
#include "operaciones_basicas.h"

double obtenerMemoriaMatrizKB(MatrixLIL* A);
double calcularResidualFinal(MatrixLIL* A, double* b, double* x);
double calcularResidualFinalCSR(MatrixCSR* A, double* b, double* x);

#endif
