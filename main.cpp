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

	int mypid = _getpid();
	//std::cout << " Rank [" << rank << "]: Pid = " << mypid << std::endl;

	ProcessController PC(PETSC_COMM_WORLD, true);
	if (rank)
	{
		PC.handle(rank, size);
		PetscFinalize();
		return 0;
	}

	// Invoking solver in helper processes
	PC.sendCommand(4);

	std::unique_ptr<SolverTest> test = std::make_unique<PETScCGTest>();
	test->argc = argc;
	test->argv = argv;
	test->A_name = "C:\\tcybakov\\MPI_Ksp_Solver\\repository\\Testing\\TestHomoStress4k.mtx";
	test->file_format = MatrixFileFormat::MTX;
	test->test();

	PC.shutdown();

	PetscFinalize();
	return 0;
}

// "C:\\tcybakov\\MPI_Ksp_Solver\\repository\\Testing\\TestHomoStress4k.mtx"
//  test->A_name = "c:\\Ulkin\\Projects\\GitHub\\MPI_T1\\workdir\\Test.mtx";

/*cout << "Number of processes: " << size << endl;
cout << "Press Enter........................" << endl;
cout.flush();
getchar();*/

//char buffer[1024];
//sprintf_s(buffer, "MyPid = |%d|\n", mypid);
//cout << buffer;

//if (rank) PetscBarrier(PETSC_NULL);