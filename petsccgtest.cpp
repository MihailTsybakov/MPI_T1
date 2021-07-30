#include "solvertest.h"

#include <process.h>
#include <petscksp.h>

using namespace std;

PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void* dummy)
{
  if (!(n % 50)) PetscPrintf(PETSC_COMM_WORLD,
    "Iteration %D KSP Residual norm %14.12e \n", n, rnorm);
  return 0;
}


int PETScCGTest::test()
{
	read_matrix_file(string(A_name), file_format, ia, ja, a, sizea, nnza);
	if (sizea == 0) return -1;

	int err = testSpecific();

	const int scatter_time = cr::duration_cast<cr::milliseconds>(scatter_end - scatter_start).count();
	const int time = cr::duration_cast<cr::milliseconds>(end - start).count();
	cout << "Scatter time == " << scatter_time << " ms" << endl; // ~~~
	cout << "Solve time   == " << time << " ms" << endl;
	cout << "Error code   == " << err << std::endl;
}

int localSize(int global_size, int mpi_rank, int mpi_size)
{
  return global_size / mpi_size + ((mpi_rank < (global_size % mpi_size)) ? 1 : 0);
}

int PETScCGTest::testSpecific()
{
	ProcessController pc(true);

	// Updating context
	pc.context.matrix_size = sizea;

	// Sending matrix size
	pc.evaluateTask(1);

	scatter_start = cr::system_clock::now();

	// Scattering IA (Everything is stored in local root context until evaluateTask() is called)
	pc.context.ia = ia; // ~~!!
	pc.context.ja = ja;
	pc.context.a = a;

	pc.context.loc_ia_sizes.resize(pc.mpiSize);
	pc.context.loc_ia_starts.resize(pc.mpiSize);

	pc.context.loc_ia_sizes[0] = localSize(sizea, 0, pc.mpiSize);
	pc.context.loc_ia_starts[0] = 1;
	for (int rank = 1; rank < pc.mpiSize; ++rank)
	{
		pc.context.loc_ia_sizes[rank] = localSize(sizea, rank, pc.mpiSize);
		pc.context.loc_ia_starts[rank] = pc.context.loc_ia_starts[rank - 1] + pc.context.loc_ia_sizes[rank - 1];
	}
	for (int rank = pc.mpiSize - 1; rank >= 0; rank--)
	{
		for (int i = 0; i < pc.context.loc_ia_sizes[rank]; ++i)
		{
			pc.context.ia[pc.context.loc_ia_starts[rank] + i] -= pc.context.ia[pc.context.loc_ia_starts[rank] - 1];
		}
	}

	pc.context.loc_rows = localSize(sizea, pc.mpiRank, pc.mpiSize);

	// Scatterv
	pc.evaluateTask(2);

	// Scattering JA
	pc.context.loc_val_sizes.resize(pc.mpiSize);  // Excess?
	pc.context.loc_val_starts.resize(pc.mpiSize); // Excess?

	int sum = 0;
	for (int rank = 0; rank < pc.mpiSize; ++rank)
	{
		// Number of elements is in the last local ia
		int last_loc_ia_idx = pc.context.loc_ia_starts[rank] + pc.context.loc_ia_sizes[rank] - 1;
		pc.context.loc_val_sizes[rank] = pc.context.ia[last_loc_ia_idx];
		pc.context.loc_val_starts[rank] = sum;
		sum += pc.context.loc_val_sizes[rank];
	}

	// Scatterv
	pc.evaluateTask(3);

	// Scattering A
	pc.evaluateTask(4);

	scatter_end = cr::system_clock::now();

	// Removing full ia, ja, a from root context:
	pc.context.ia.resize(0);
	pc.context.ja.resize(0);
	pc.context.a.resize(0);

	// Converting to Petsc data types:
	pc.evaluateTask(7);

	// Creating MPI matrix
	pc.evaluateTask(8);

	// Releasing Petsc memory
	pc.evaluateTask(9);

	// MPI matrix assembly
	pc.evaluateTask(10);

	// Block A: Forming vectors
	pc.evaluateTask(11);

	// Block B: Solver set up
	pc.evaluateTask(12);

	start = cr::system_clock::now();

	// Invoking solver
	pc.evaluateTask(13);

	end = cr::system_clock::now();

	// Checking solution params & convergence
	pc.evaluateTask(14);

	cout << "The system has been solved\n";
	cout << "Number of iterations: " << pc.context.iterations << endl;
	cout << "Residual norm: " << pc.context.res_norm << endl;

	// Releasing solver memory
	pc.evaluateTask(15);

	// Shutting down all processes:
	pc.evaluateTask(-1);
	
	return 0;  
}
