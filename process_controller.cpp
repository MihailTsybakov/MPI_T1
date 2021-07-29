#include "process_controller.h"

ProcessController::ProcessController(MPI_Comm communicator, bool logs) : communicator(communicator), logs(logs)  {}

void ProcessController::handle(int rank, int size)
{
	bool handle_flag = true;
	while (handle_flag)
	{
		int command;
		MPI_Bcast(&command, 1, MPI_INT, 0, communicator);

		switch (command)
		{
		case 0:
			// Finalizing handle
			if (logs) std::cout << " <logs> Rank [" << rank << "]: Shutting down. " << std::endl;
			handle_flag = false;
			break;
		case 1:
			// Pi
			double PI;
			calc_PI(1e7, 1e-7, rank, size, PI);
			break;
		case 2:
			// Area
			double circle_area;
			calc_area(1e7, 1e-7, 3.0, rank, size, circle_area);
			break;
		case 3:
			// Volume
			double sphere_volume;
			calc_volume(1e7, 1e-7, 3.0, rank, size, sphere_volume);
			break;
		case 4:
		{
			// Solver call
			std::unique_ptr<SolverTest> test = std::make_unique<PETScCGTest>();
			test->argc = 0; //~~
			test->argv = nullptr; //!!
			test->A_name = "C:\\tcybakov\\MPI_Ksp_Solver\\repository\\Testing\\TestHomoStress4k.mtx";
			test->file_format = MatrixFileFormat::MTX;
			test->test();
			break;
		}
		default:
			// Unknown command
			std::cout << " Rank [" << rank << "]: Unknown command (" << command << ") recieved. Shutting down." << std::endl;
			handle_flag = false;
			break;
		}
	}
	//MPI_Finalize();
}

void ProcessController::sendCommand(int command)
{
	MPI_Bcast(&command, 1, MPI_INT, 0, communicator);
}

void ProcessController::shutdown()
{
	sendCommand(0);
}