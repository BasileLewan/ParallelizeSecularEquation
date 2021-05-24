OBJS	= main.o util.o benchmark.o main_omp.o util_omp.o benchmark_omp.o
OUT	= main benchmark

OBJS0	= main.o util.o
SOURCE0	= main.c util.c
HEADER0	= util.h
OUT0	= main

OBJS1	= benchmark.o main.o util.o
SOURCE1	= benchmark.c main.c util.c
HEADER1	= main.h util.h
OUT1	= benchmark

OBJS2	= benchmark_omp.o main_omp.o util_omp.o
SOURCE2	= benchmark.c main.c util.c
HEADER2	= main.h util.h
OUT2	= benchmark

OBJS3	= main_omp.o util_omp.o
SOURCE3	= main.c util.c
HEADER3	= main.h util.h
OUT3	= main

CC	 = gcc-10
FLAGS	 = -g -c -Wall -Wextra -O3 -ftree-vectorize
LFLAGS	 = -lm

all: main_omp benchmark_omp

main: $(OBJS0) $(LFLAGS)
	$(CC) -g $(OBJS0) -o $(OUT0)

benchmark: $(OBJS1) $(LFLAGS)
	$(CC) -g $(OBJS1) -o $(OUT1)
	
benchmark_omp: $(OBJS2) $(LFLAGS)
	$(CC) -g $(OBJS2) -lm -fopenmp -o $(OUT2)

main_omp: $(OBJS3) $(LFLAGS)
	$(CC) -g $(OBJS3) -lm -fopenmp -o $(OUT3)

main.o: main.c
	$(CC) $(FLAGS) main.c 

util.o: util.c,
	$(CC) $(FLAGS) util.c

benchmark.o: benchmark.c
	$(CC) $(FLAGS) bechmark.c
	
main_omp.o: main.c
	$(CC) $(FLAGS) -fopenmp main.c -o main_omp.o

util_omp.o: util.c
	$(CC) $(FLAGS) -fopenmp util.c -o util_omp.o

benchmark_omp.o: benchmark.c
	$(CC) $(FLAGS) -fopenmp benchmark.c -o benchmark_omp.o



clean:
	rm -f $(OBJS) $(OUT)