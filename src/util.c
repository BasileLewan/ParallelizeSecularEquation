#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "util.h"
#define SIZE 100

const extern size_t N = SIZE;

double** add(double** m1, double** m2) {
	double** res = alloc2d(N);

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++)
			res[i][j] = m1[i][j] + m2[i][j];
	}
	return res;
}

double sqr(double x) {
	return x * x;
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
			m[0][k] *= x; // Violation d'acces lors de la lecture
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
		printf("�|\n");
	}
	for (unsigned int k = 0; k < N; k++)
		printf("----------");
	printf("\n");
}

void print1(double* m) {	//print 1st dimension arrays
	printf("{");
	for(int i = 0; i < N - 1; ++i) {
		printf("%03.1f, ", m[i]);
	}
	printf("%03.1f}\n", m[N - 1]);
}

double** gen_sym() {
	/* Generate random symetric matrix */
	double** res = alloc2d(N);
	double a;
	srand(time(NULL));

	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = 0; j < i + 1; j++) {
			a = rand() % 10000;
			res[i][j] = a;
			res[j][i] = a;
		}
	}
	return res;
}

double** eye() {
	double** res = alloc2d(N);
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++)
			res[i][j] = i == j ? 1 : 0;
	}
	return res;
}

double** tridiag(){
	/* Householder's tridiagonalization */
	unsigned i;
	double s = 0, a;
	double* v = malloc(N * sizeof(double));
	for (int i = 0; i < N; ++i)
		v[i] = 0;
	double** I = eye();
	double** A_old = gen_sym();
	double** tmp; double** P; double** A_new = NULL;

	for (int k = 0; k < N - 2; k++) {
		//print(A_old);
		for (i = 0; i < k + 1; ++i)
			v[i] = 0;
		for (i = k + 1; i < N; i++)
			s += sqr(A_old[i][k]);
		s = sqrt(s);
		v[k + 1] = sqrt(.5 * (1 + fabs(A_old[k + 1][k]) / (s + 2 * DBL_EPSILON)));
		double sgn = (double)(A_old[k + 1][k] > 0) - (double)(A_old[k + 1][k] < 0);
		for (i = k + 2; i < N; i++)
			v[i] = A_old[i][k] * sgn / (2 * v[k + 1] * s + 2* DBL_EPSILON);
		tmp = alloc2d(N);
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < i + 1; j++) {
				a = -2 * v[i] * v[j];
				tmp[i][j] = a;
				tmp[j][i] = a;
			}
		}
		P = add(I, tmp);
		free(tmp[0]);
		tmp = multiply(P, A_old);
		A_new = multiply(tmp, P);
		// print(P); print(A_new);
		free(tmp[0]); free(P[0]); free(A_old[0]);
		A_old = A_new;
	}
	return A_new;
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
