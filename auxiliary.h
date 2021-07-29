#ifndef AUX
#define AUX

#include <iostream>
#include <thread>
#include <chrono>

#include "mpi.h"

/// Calculates PI
void calc_PI(int N, double d, int rank, int size, double& PI);

/// Calculates circle area
void calc_area(int N, double d, double R, int rank, int size, double& circle_area);

/// Calculates sphere volume
void calc_volume(int N, double d, double R, int rank, int size, double& sphere_vol);

//int old_main(int argc, char* argv[]);

#endif//AUX
