#include <stdlib.h>
#include <stdio.h>
#include "util.h"
#define SIZE 4

const size_t N = SIZE;

double** add(double** m1, double** m2) {
	double** res = alloc2d(N);

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) 
			res[i][j] = m1[i][j] + m2[i][j];	
	}
	return res;
}

double** multiply(double** m1, double** m2)
{
	double** res = alloc2d(N);
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			double tmp = 0;
			int k;
#pragma omp parallel for reduction(+:tmp) 
			for (k = 0; k < N; k++) {
				tmp += m1[i][k] * m2[k][j];
			}
			res[i][j] = tmp;
		}
	}

	return res;
}

void times(double x, double** m) {
	for (unsigned int k = 0; k < N * N; k++) {
			m[0][k] *= x; // Violation d'accès lors de la lecture 
	}
}

void print(double** m) {
	for (unsigned int k = 0; k < N; k++)
		printf("----------");
	printf("\n");
	for (unsigned int i = 0; i < N; i++) {
		printf("|  ");
		for (unsigned int j = 0; j < N; j++)
			printf("%07.1f  ", m[i][j]);
		printf(" |\n");
	}
	for (unsigned int k = 0; k < N; k++)
		printf("----------");
	printf("\n");
}

double** rd_gen() {
	double** res = alloc2d(N);
	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = 0; j < N; j++)
			res[i][j] = rand() % 10000;
	}
	return res;
}

double** alloc2d(unsigned sz) {
	double** res = (double**)malloc(sz * sizeof(double*)); // Exception ici ???
	if (res == NULL)
		perror("Memory error");
	if ((res[0] = malloc(sz * sz * sizeof(double))) == NULL)
		perror("Memory error");
	for (unsigned i = 1; i < sz; i++)
		res[i] = res[0] + i * sz;

	return (double**)res;
}