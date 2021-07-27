#include <iostream>
#include "mpi.h"

int main(int argc, char* argv[])
{
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int REPEATS = 1e7;
	double d = 1e-7;
	double d2 = 1e-14;

	double I = 0.0;
	double result = 0.0;

	for (int i = rank; i < REPEATS; i += size)
	{
		I += 1 / (1 + d2 * i * i);
	}

	MPI_Reduce(&I, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (!rank) std::cout << "Pi = " << result * 4 * d << std::endl;

	MPI_Finalize();

	return 0;
}