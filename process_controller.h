#ifndef PROC_CONTR
#define PROC_CONTR

#include "mpi.h"
#include "auxiliary.h"
#include "solvertest.h"
#include <petscksp.h>

class ProcessTask
{
private:

public:
	ProcessTask() {};
	virtual void task() = 0;
	virtual ~ProcessTask() {};
};



class ProcessController
{
private:
	std::map<int, std::function<void(...)>> tasks;
	MPI_Comm communicator = PETSC_COMM_WORLD;
	bool logs;
public:
	ProcessController(bool logs = false);
	void handle_tmp(int command = -1);
	void handle(ProcessTask& task, int command = -1);
	void sendCommand(int command);
	void shutdown();
	void getInfo(int& rank, int& size) const;
	void wait();
	void yield();
};

#endif//PROC_CONTR
