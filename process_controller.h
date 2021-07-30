#ifndef PROC_CONTR
#define PROC_CONTR

#include "mpi.h"
#include "auxiliary.h"
#include "solvertest.h"

#include <petscksp.h>
#include <string>

class ProcessTask
{
private:

public:
	ProcessTask() {};
	virtual void task() = 0;
	virtual ~ProcessTask() {};
};

//class ProcessController;
class ProcessContext
{
public:
	// Common for all processes:
	int matrix_size;
	std::vector<int> loc_ia_sizes;
	std::vector<int> loc_ia_starts;
	std::vector<int> loc_val_sizes;
	std::vector<int> loc_val_starts;

	// Not synchronized, local context :
	int loc_rows;
	int loc_num;
	Mat A;
	Vec result, ref_result, b;
	KSP ksp;
	PC pc;

	PetscInt iterations;
	PetscScalar res_norm;

	std::vector<int> loc_ja;
	std::vector<int> loc_ia;  
	std::vector<double> loc_a;
	
	PetscInt* petsc_loc_ia;
	PetscInt* petsc_loc_ja;
	PetscScalar* petsc_loc_a;

	std::vector<int> ia; // For root only (may be readed here)
	std::vector<int> ja; // For root only (may be readed here)
	std::vector<double> a; // For root only (may be readed here)
	
	//friend class ProcessController;
};

class ProcessController
{
private:
	std::map<int, std::function<void(double)>> tasks;
	MPI_Comm communicator = PETSC_COMM_WORLD;
	bool logs;
public:
	ProcessContext context;
	int mpiRank, mpiSize;

	ProcessController(bool logs = false);
	void handle_tmp(int command = -1);
	void handle(ProcessTask& task, int command = -1);
	void sendCommand(int command);
	void shutdown();
	void getInfo(int& rank, int& size) const;
	void wait();
	void yield();

	void syncContext();
	void evaluateTask(int taskID);
	void waitForTask();
};


int localMatrixSize(int matrix_size, int mpi_process_rank, int mpi_group_size);
PetscErrorCode MyKSPMonitor_(KSP ksp, PetscInt n, PetscReal rnorm, void* dummy);
#endif//PROC_CONTR
