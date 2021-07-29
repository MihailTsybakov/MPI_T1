#ifndef PROC_CONTR
#define PROC_CONTR

#include "mpi.h"
#include "auxiliary.h"
#include "solvertest.h"

class ProcessController
{
private:
	MPI_Comm communicator;
	bool logs;
public:
	ProcessController(MPI_Comm communicator, bool logs);
	void handle(int rank, int size);
	void sendCommand(int command);
	void shutdown();
};

#endif//PROC_CONTR
