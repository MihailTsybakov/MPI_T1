#include "auxiliary.h"

/// Calculates PI
void calc_PI(int N, double d, int rank, int size, double& PI)
{
	std::cout << " > rank [" << rank << "] working on PI" << std::endl;
	double I = 0.0, res = 0.0;
	for (int i = rank; i < N; i += size)
	{
		I += d / (1 + d * d * i * i);
	}
	MPI_Reduce(&I, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	PI = 4 * res;
}

/// Calculates circle area
void calc_area(int N, double d, double R, int rank, int size, double& circle_area)
{
	std::cout << " > rank [" << rank << "] working on AREA" << std::endl;
	double I = 0.0, res = 0.0;
	for (int i = rank; i < N; i += size)
	{
		double h = i*d*R;
		I += std::sqrt(R * R - h * h) * 2;
	}
	MPI_Reduce(&I, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	circle_area = 2 * R*d*res;
}

/// Calculates sphere volume
void calc_volume(int N, double d, double R, int rank, int size, double& sphere_vol)
{
	std::cout << " > rank [" << rank << "] working on VOLUME" << std::endl;
	double I = 0.0, res = 0.0;
	for (int i = rank; i < N; i += size)
	{
		double h = i * d * R;
		I += 3.1415*(R * R - h * h);
	}
	MPI_Reduce(&I, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	sphere_vol = 2 * R*d*res;
}

//int old_main(int argc, char* argv[])
//{
//	int rank, size;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//	process_controller PC(MPI_COMM_WORLD);
//
//	if (rank)
//	{
//		PC.handle(rank, size);
//		return 0;
//	}
//
//	// Doing something in root...
//	std::cout << std::endl;
//
//	std::cout << "Root alone 1.1" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	std::cout << "Root alone 1.2" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	std::cout << "Root alone 1.3" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//	// Done doing something, calculating PI now
//
//	double PI;
//	int root_command = 1;
//
//	MPI_Bcast(&root_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	calc_PI(1e7, 1e-7, rank, size, PI);
//
//	// Doing something in root again
//	std::cout << "Calculated PI: " << PI << std::endl;
//
//	std::cout << "Root alone 2.1" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	std::cout << "Root alone 2.2" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	std::cout << "Root alone 2.3" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	// Done doing something, calculating circle area:
//
//	double circle_area;
//	root_command = 2;
//
//	MPI_Bcast(&root_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	calc_area(1e7, 1e-7, 3.0, rank, size, circle_area);
//
//	// Doing something in root again
//	std::cout << "Calaculated area: " << circle_area << std::endl;
//
//	std::cout << "Root alone 3.1" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	std::cout << "Root alone 3.2" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//
//	std::cout << "Root alone 3.3" << std::endl;
//	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
//	// Done doing something, calculating sphere volume
//
//	double sphere_vol;
//	root_command = 3;
//
//	MPI_Bcast(&root_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	calc_volume(1e7, 1e-7, 3.0, rank, size, sphere_vol);
//
//	std::cout << "Calculated volume: " << sphere_vol << std::endl;
//
//	// Shutting down
//
//	root_command = 0;
//	MPI_Bcast(&root_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//	MPI_Finalize();
//	return 0;
//}