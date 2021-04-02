#include <stdlib.h>


extern const size_t N;

/* Basical matrix operations*/
double** rd_gen();
double** add(double** m1, double** m2);
double** multiply(double** m1, double** m2);
void times(double x, double** m);
void print(double** m);
double** alloc2d(unsigned size);