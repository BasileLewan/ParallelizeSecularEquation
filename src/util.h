#pragma once
#include <stdlib.h>


extern size_t N;
double sqr(double x);

/* Basical matrix operations*/
double** add(double** m1, double** m2);
double** multiply(double** m1, double** m2);
void times(double x, double** m);
void print(double** m);

/* Matrices generation */
double** alloc2d(unsigned size);
double** gen_sym();
double** eye();

/* Serious business */
double** tridiag(size_t size);
double get_rho(double** matrix);
double* get_ksi(double** matrix);
double* get_delta(double** matrix, double* ksi, double rho);