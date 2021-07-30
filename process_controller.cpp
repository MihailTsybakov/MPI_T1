#include "process_controller.h"

ProcessController::ProcessController(bool logs) : logs(logs)
{
	getInfo(mpiRank, mpiSize);
	context.loc_ia_sizes.resize(mpiSize); //~~
	context.loc_ia_starts.resize(mpiSize);
	context.loc_val_sizes.resize(mpiSize);
	context.loc_val_starts.resize(mpiSize);//~~
	// -1 = Shutdown function
	tasks[-1] = [&](double tmp)
	{
		if (logs) std::cout << " <logs> Rank [" << mpiRank << "]: Shutting down." << std::endl;
		PetscFinalize();
	};
	// 0 = Sync Context
	tasks[0] = [&](double tmp)
	{
		syncContext();
	};
	// 1 = Send matrix size
	tasks[1] = [this](double tmp)
	{
		MPI_Bcast(&(context.matrix_size), 1, MPI_INT, 0, communicator);
	};
	// 2 = ScatterV_IA
	tasks[2] = [this](double tmp)
	{
		context.loc_rows = localMatrixSize(context.matrix_size, mpiRank, mpiSize);
		context.loc_ia.resize(context.loc_rows + 1);
		context.loc_ia[0] = 0;
		// ~~!!
		MPI_Scatterv(context.ia.data(), context.loc_ia_sizes.data(), context.loc_ia_starts.data(), MPI_INT,
			context.loc_ia.data() + 1, context.loc_rows, MPI_INT, 0, communicator);
	};
	// 3 = ScatterV_JA
	tasks[3] = [this](double tmp)
	{
		context.loc_num = context.loc_ia[context.loc_ia.size() - 1];
		context.loc_ja.resize(context.loc_num);
		MPI_Scatterv(context.ja.data(), context.loc_val_sizes.data(), context.loc_val_starts.data(), MPI_INT,
			context.loc_ja.data(), context.loc_num, MPI_INT, 0, communicator);
	};
	// 4 = ScatterV_A
	tasks[4] = [this](double tmp)
	{
		context.loc_a.resize(context.loc_num);
		MPI_Scatterv(context.a.data(), context.loc_val_sizes.data(), context.loc_val_starts.data(), MPI_DOUBLE,
			context.loc_a.data(), context.loc_num, MPI_DOUBLE, 0, communicator);
	};
	// 5 = Resize local info vectors
	tasks[5] = [this](double tmp)
	{
		if (mpiRank)
		{
			context.loc_ia_sizes.resize(mpiSize);
			context.loc_ia_starts.resize(mpiSize);
			context.loc_val_sizes.resize(mpiSize);
			context.loc_val_starts.resize(mpiSize);
		}
	};
	// 6 = Debugging task: prints context's content for each process
	tasks[6] = [this](double tmp)
	{
		std::cout << " Rank [" << mpiRank << "]: -------------------------" << std::endl;

		std::cout << " Rank [" << mpiRank << "]: Matrix size = " << context.matrix_size << std::endl;
		std::cout << " Rank [" << mpiRank << "]: Local rows = " << context.loc_rows << std::endl;
		std::cout << " Rank [" << mpiRank << "]: Local ia sizes: ";
		for (int i = 0; i < mpiSize; ++i) std::cout << context.loc_ia_sizes[i] << " ";
		std::cout << std::endl;
		std::cout << " Rank [" << mpiRank << "]: My ja local part: ";
		for (int i = 0; i < context.loc_num; ++i) std::cout << context.loc_ja[i] << " ";
		std::cout << std::endl;
		std::cout << " Rank [" << mpiRank << "]: My ia local part: ";
		for (int i = 0; i < context.loc_rows+1; ++i) std::cout << context.loc_ia[i] << " ";
		std::cout << std::endl;
		std::cout << " Rank [" << mpiRank << "]: -------------------------" << std::endl;
	};
	// 7 = Allocating memory for Petsc types and copying data
	tasks[7] = [this](double tmp)
	{
		PetscMalloc(context.loc_ia.size() * sizeof(PetscInt), &context.petsc_loc_ia);
		PetscMalloc(context.loc_ja.size() * sizeof(PetscInt), &context.petsc_loc_ja);
		PetscMalloc(context.loc_a.size() * sizeof(PetscScalar), &context.petsc_loc_a);
		for (size_t i = 0; i < context.loc_ia.size(); ++i) context.petsc_loc_ia[i] = context.loc_ia[i];
		for (size_t i = 0; i < context.loc_ja.size(); ++i) context.petsc_loc_ja[i] = context.loc_ja[i];
		for (size_t i = 0; i < context.loc_a.size(); ++i) context.petsc_loc_a[i] = context.loc_a[i];
	};
	// 8 = Matrix create
	tasks[8] = [this](double tmp)
	{
		MatCreateMPIAIJWithArrays(communicator, context.loc_rows, PETSC_DECIDE,
			PETSC_DETERMINE, context.matrix_size, context.petsc_loc_ia, context.petsc_loc_ja,
			context.petsc_loc_a, &context.A);
	};
	// 9 = Petsc memory releasing
	tasks[9] = [this](double tmp)
	{
		PetscFree(context.petsc_loc_a);
		PetscFree(context.petsc_loc_ia);
		PetscFree(context.petsc_loc_ja);
	};
	// 10 = Matrix assembly
	tasks[10] = [this](double tmp)
	{
		MatAssemblyBegin(context.A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(context.A, MAT_FINAL_ASSEMBLY);
	};
	// 11 = Forming vectors
	tasks[11] = [this](double tmp)
	{
		MatCreateVecs(context.A, &context.b, NULL);
		VecDuplicate(context.b, &context.result); 
		VecDuplicate(context.b, &context.ref_result);
		VecSet(context.result, 0.0);
		VecSet(context.ref_result, 1.0);
		MatMult(context.A, context.ref_result, context.b);
	};
	// 12 = Setting up solver
	tasks[12] = [this](double tmp)
	{
		KSPCreate(communicator, &context.ksp);
		KSPGetPC(context.ksp, &context.pc);
		PCSetType(context.pc, PCJACOBI);
		KSPMonitorSet(context.ksp, MyKSPMonitor_, NULL, 0);
		KSPSetOperators(context.ksp, context.A, context.A);
		//KSPSetType(ksp, KSPCG);
		KSPSetType(context.ksp, KSPGMRES);
		KSPSetTolerances(context.ksp, 1e-8, PETSC_DEFAULT,
			PETSC_DEFAULT, 1'000);
		KSPSetFromOptions(context.ksp);

		// These calls are optional enable more precise profiling, 
		// since both will be called within KSPSolve() if they 
		// haven't been called already.
		KSPSetUp(context.ksp);
	};
	// 13 = Solver invoke
	tasks[13] = [this](double tmp)
	{
		KSPSolve(context.ksp, context.b, context.result);
	};
	// 14 = Convergence check
	tasks[14] = [this](double tmp)
	{
		KSPConvergedReason reason;
		KSPGetConvergedReason(context.ksp, &reason);
		if (reason < 0) {
			cout << "Divergence detected: " << reason << endl;
		}
		else
		{
			PetscInt its;
			KSPGetIterationNumber(context.ksp, &its);
			PetscScalar res_norm;
			KSPGetResidualNorm(context.ksp, &res_norm);
			context.iterations = its;
			context.res_norm = res_norm;
		}
	};
	// 15 = Solver memory releasing
	tasks[15] = [this](double tmp)
	{
		KSPDestroy(&context.ksp);
		MatDestroy(&context.A);
		VecDestroy(&context.ref_result);
		VecDestroy(&context.result);
		VecDestroy(&context.b);
	};
}

void ProcessController::syncContext()
{
	MPI_Bcast(&(context.matrix_size), 1, MPI_INT, 0, communicator);

	MPI_Bcast(&(context.loc_ia_sizes.at(0)), mpiSize, MPI_INT, 0, communicator);
	MPI_Bcast(&(context.loc_ia_starts.at(0)), mpiSize, MPI_INT, 0, communicator);
	MPI_Bcast(&(context.loc_val_sizes.at(0)), mpiSize, MPI_INT, 0, communicator);
	MPI_Bcast(&(context.loc_val_starts.at(0)), mpiSize, MPI_INT, 0, communicator);

}

// Specially for root process
void ProcessController::evaluateTask(int taskID)
{
	double service_arg;
	if (!mpiRank)
	{
		syncContext();
		MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);
	}
	tasks[taskID](service_arg);
	if (logs) std::cout << " <logs> Rank [" << mpiRank << "]: Task " << taskID << " done." << std::endl;
}

// Specially for helper processes
void ProcessController::waitForTask()
{
	bool finish_flag = false;
	while (!finish_flag)
	{
		int taskID;
		syncContext();
		MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);
		if (taskID == -1)
		{
			finish_flag = true;
		}
		evaluateTask(taskID);
	}
}

void ProcessController::handle_tmp(int command)
{
	bool handle_flag = true;
	while (handle_flag)
	{
		if (mpiRank)
		{
			MPI_Bcast(&command, 1, MPI_INT, 0, communicator);
		}
		else
		{
			sendCommand(command);
			handle_flag = false;
		}

		//tasks[command];

		switch (command)
		{
		case 0:
			// Finalizing handle
			if (logs) std::cout << " <logs> Rank [" << mpiRank << "]: Shutting down. " << std::endl;
			handle_flag = false;
			PetscFinalize();
			break;
		case 1:
			// Pi
			double PI;
			calc_PI(1e7, 1e-7, mpiRank, mpiSize, PI);
			break;
		case 2:
			// Area
			double circle_area;
			calc_area(1e7, 1e-7, 3.0, mpiRank, mpiSize, circle_area);
			break;
		case 3:
			// Volume
			double sphere_volume;
			calc_volume(1e7, 1e-7, 3.0, mpiRank, mpiSize, sphere_volume);
			break;
		case 4:
		{
			// Solver call
			std::cout << " Rank [" << mpiRank << "]: Called solver" << std::endl;
			break;
		}
		case 42:

		default:
			// Unknown command
			std::cout << " Rank [" << mpiRank << "]: Unknown command (" << command << ") recieved. Shutting down." << std::endl;
			handle_flag = false;
			break;
		}
	}
}

void ProcessController::wait()
{
	MPI_Barrier(communicator);
}

void ProcessController::yield()
{
	MPI_Barrier(communicator);
}

void ProcessController::handle(ProcessTask& task, int command)
{
	int rank, size;
	getInfo(rank, size);
	bool work_flag = true;
	while (work_flag)
	{
		if (rank)
		{
			MPI_Bcast(&command, 1, MPI_INT, 0, communicator);
		}
		else
		{
			sendCommand(command);
			work_flag = false;
		}
		task.task();
	}
}


void ProcessController::sendCommand(int command)
{
	MPI_Bcast(&command, 1, MPI_INT, 0, communicator);
}

void ProcessController::shutdown()
{
	sendCommand(0);
}

void ProcessController::getInfo(int& rank, int& size) const
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
}

int localMatrixSize(int matrix_size, int mpi_process_rank, int mpi_group_size)
{
	return (matrix_size / mpi_group_size + ((mpi_process_rank < (matrix_size % mpi_group_size)) ? 1 : 0));
}

PetscErrorCode MyKSPMonitor_(KSP ksp, PetscInt n, PetscReal rnorm, void* dummy)
{
	if (!(n % 50)) PetscPrintf(PETSC_COMM_WORLD,
		"Iteration %D KSP Residual norm %14.12e \n", n, rnorm);
	return 0;
}