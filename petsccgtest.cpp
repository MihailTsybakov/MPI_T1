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

	class CGTask : public ProcessTask
	{
	public:
		void task() override
		{
			ProcessController pc;
			int rank, size;
			pc.getInfo(rank, size);
			std::cout << " >>>>>>>>>>>>>> Rank [" << rank << "]: CGTask called <<<<<<<<<<<<<" << std::endl;
		}
	};

	CGTask SomeTask;

	pc.handle(SomeTask);

  //scatter_start = cr::system_clock::now();
  /*
  int root = 0;
  // Sending matrix size here
  MPI_Bcast(&sizea, 1, MPI_INT, root, PETSC_COMM_WORLD);
 
  // make local CSR from global CSR
  // https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateMPIAIJWithArrays.html
  // ia on each processors starts with zero
  vector<int> loc_ia_sizes(mpi_size);
  vector<int> loc_ia_starts(mpi_size);
  if (!mpi_rank)
  {
    loc_ia_sizes[0] = localSize(sizea, 0, mpi_size);
    loc_ia_starts[0] = 1;
    for (int rank = 1; rank < mpi_size; ++rank)
    {
      loc_ia_sizes[rank] = localSize(sizea, rank, mpi_size);
      loc_ia_starts[rank] = loc_ia_starts[rank - 1] + loc_ia_sizes[rank - 1];
    }
    for (int rank = mpi_size - 1; rank >= 0; rank--)
    {
      for (int i = 0; i < loc_ia_sizes[rank]; ++i)
      {
        ia[loc_ia_starts[rank] + i] -= ia[loc_ia_starts[rank] - 1];
      }
    }
  }
  int loc_rows = localSize(sizea, mpi_rank, mpi_size);
  vector<int> loc_ia(loc_rows + 1);
  loc_ia[0] = 0;
  MPI_Scatterv(ia.data(), loc_ia_sizes.data(), loc_ia_starts.data(), MPI_INT,
    loc_ia.data() + 1, loc_rows, MPI_INT, root, PETSC_COMM_WORLD);

  // ja, a
  vector<int> loc_val_sizes(mpi_size);
  vector<int> loc_val_starts(mpi_size);
  if (!mpi_rank)
  {
    int sum = 0;
    for (int rank = 0; rank < mpi_size; ++rank)
    {
      // number of elements is in the last local ia
      int last_loc_ia_idx = loc_ia_starts[rank] + loc_ia_sizes[rank] - 1;
      loc_val_sizes[rank] = ia[last_loc_ia_idx];
      loc_val_starts[rank] = sum;
      sum += loc_val_sizes[rank];
    }
  }
  int loc_num = loc_ia[loc_ia.size() - 1];
  vector<int> loc_ja(loc_num);
  MPI_Scatterv(ja.data(), loc_val_sizes.data(), loc_val_starts.data(), MPI_INT,
    loc_ja.data(), loc_num, MPI_INT, root, PETSC_COMM_WORLD);
  vector<double> loc_a(loc_num);
  MPI_Scatterv(a.data(), loc_val_sizes.data(), loc_val_starts.data(), MPI_DOUBLE,
    loc_a.data(), loc_num, MPI_DOUBLE, root, PETSC_COMM_WORLD);

  // remove global
  if (!mpi_rank) ia.resize(0); ja.resize(0); a.resize(0);

  PetscBarrier(PETSC_NULL);
  if (!mpi_rank) scatter_end = cr::system_clock::now();

  // convert int vector to bigger memory size PetscInt array
  PetscInt* petsc_loc_ia;
  PetscInt* petsc_loc_ja;
  PetscScalar* petsc_loc_a;
  ierr = PetscMalloc(loc_ia.size() * sizeof(PetscInt), &petsc_loc_ia); CHKERRQ(ierr);
  ierr = PetscMalloc(loc_ja.size() * sizeof(PetscInt), &petsc_loc_ja); CHKERRQ(ierr);
  ierr = PetscMalloc(loc_a.size() * sizeof(PetscScalar), &petsc_loc_a); CHKERRQ(ierr);
  for (int i = 0; i < loc_ia.size(); ++i) petsc_loc_ia[i] = static_cast<PetscInt>(loc_ia[i]);
  for (int i = 0; i < loc_ja.size(); ++i) petsc_loc_ja[i] = static_cast<PetscInt>(loc_ja[i]);
  for (int i = 0; i < loc_a.size(); ++i) petsc_loc_a[i] = static_cast<PetscInt>(loc_a[i]);

  Mat A;
  ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, loc_rows, PETSC_DECIDE,
    PETSC_DETERMINE, sizea, petsc_loc_ia, petsc_loc_ja,
    petsc_loc_a, &A); CHKERRQ(ierr);

  ierr = PetscFree(petsc_loc_ia); CHKERRQ(ierr);
  ierr = PetscFree(petsc_loc_ja); CHKERRQ(ierr);
  ierr = PetscFree(petsc_loc_a); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec result, ref_result, b;
  ierr = MatCreateVecs(A, &b, NULL); CHKERRQ(ierr);
  ierr = VecDuplicate(b, &result); CHKERRQ(ierr);
  ierr = VecDuplicate(b, &ref_result); CHKERRQ(ierr);
  ierr = VecSet(result, 0.0); CHKERRQ(ierr);
  ierr = VecSet(ref_result, 1.0); CHKERRQ(ierr);
  ierr = MatMult(A, ref_result, b); CHKERRQ(ierr);
//  ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  KSP ksp;
  PC pc;
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
  ierr = KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
//  ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, tol, PETSC_DEFAULT,
    PETSC_DEFAULT, maxiter); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  // These calls are optional enable more precise profiling, 
  // since both will be called within KSPSolve() if they 
  // haven't been called already.
  ierr = KSPSetUp(ksp); CHKERRQ(ierr);

  PetscBarrier((PetscObject)A);
  if (!mpi_rank) start = cr::system_clock::now();
//  VecView(b, PETSC_VIEWER_STDOUT_WORLD);
//  MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  ierr = KSPSolve(ksp, b, result); CHKERRQ(ierr);

  PetscBarrier((PetscObject)A);
  if (!mpi_rank) end = cr::system_clock::now();

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  if (reason < 0) {
    cout << "Divergence detected: " << reason << endl;
  }
  else
  {
    PetscInt its;
    ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
    PetscScalar res_norm;
    ierr = KSPGetResidualNorm(ksp, &res_norm); CHKERRQ(ierr);
    if (!mpi_rank)
    {
      cout << "The system has been solved\n";
      cout << "Number of iterations: " << its << endl;
      cout << "Residual norm: " << res_norm << endl;
    }
  }

  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&ref_result); CHKERRQ(ierr);
  ierr = VecDestroy(&result); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  */
	
	
	return 0;  
}
