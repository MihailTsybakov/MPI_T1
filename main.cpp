#include <iostream>
#include <thread>
#include <chrono>

#include "auxiliary.h"
#include "process_controller.h"
#include "solvertest.h"

#include <petscksp.h>
#include <process.h>

int main(int argc, char* argv[])
{
	PetscErrorCode ierr;
	PetscMPIInt rank, size;

	ierr = PetscInitialize(&argc, &argv, (char*)0, (char*)0);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

#ifdef _DEBUG
	int mypid = _getpid();
	std::cout << " Rank [" << rank << "]: Pid = " << mypid << std::endl;
	if (!rank) 
	{
		std::cout << "Press any key to continue..." << std::endl;
		getchar();
	}
	PetscBarrier(PETSC_NULL);
#endif

	ProcessController pc(true);
	if (rank)
	{
		pc.waitForTask();
		return 0;
	}
	

	std::unique_ptr<SolverTest> test = std::make_unique<PETScCGTest>();
	test->argc = argc;
	test->argv = argv;
	test->A_name = "C:\\tcybakov\\MPI_Ksp_Solver\\repository\\Testing\\Test.mtx";
	test->file_format = MatrixFileFormat::MTX;
	test->test();

	//PetscFinalize(); //~~
	//pc.handle_tmp(0);
	
	return 0;
}

// "C:\\tcybakov\\MPI_Ksp_Solver\\repository\\Testing\\TestHomoStress4k.mtx"