#include "process_controller.h"

ProcessController::ProcessController(bool logs) : logs(logs)  {}

void ProcessController::handle_tmp(int command)
{
	int rank, size;
	getInfo(rank, size);
	bool handle_flag = true;
	while (handle_flag)
	{
		if (rank)
		{
			MPI_Bcast(&command, 1, MPI_INT, 0, communicator);
		}
		else
		{
			sendCommand(command);
			handle_flag = false;
		}

		tasks[command];

		switch (command)
		{
		case 0:
			// Finalizing handle
			if (logs) std::cout << " <logs> Rank [" << rank << "]: Shutting down. " << std::endl;
			handle_flag = false;
			PetscFinalize();
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
			std::cout << " Rank [" << rank << "]: Called solver" << std::endl;
			break;
		}
		default:
			// Unknown command
			std::cout << " Rank [" << rank << "]: Unknown command (" << command << ") recieved. Shutting down." << std::endl;
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