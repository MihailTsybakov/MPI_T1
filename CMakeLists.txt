cmake_minimum_required(VERSION 3.1)

project(MPI_Integral)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

set(SOLVTEST_H 
	solvertest.h
	auxiliary.h
	process_controller.h
)

set(SOLVTEST_SRC 
	main.cpp
	solvertest.cpp
	petsccgtest.cpp
	auxiliary.cpp
	process_controller.cpp
)

find_package(MPILibraries REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDES})
link_directories(${MKL_LIBRARYDIR})
link_directories(${MKL_BINARYDIR})
link_directories(${OMP_LIBRARYDIR})

find_package(PETSc REQUIRED)

add_subdirectory(matrixio)

add_executable(${PROJECT_NAME} ${SOLVTEST_H} ${SOLVTEST_SRC})
target_include_directories(${PROJECT_NAME} PRIVATE ${PETSC_INCLUDES} ${PETSC_CONF_INCLUDES_64})
target_link_libraries(${PROJECT_NAME} LINK_PRIVATE ${PETSC_LIBRARY_64})
target_link_libraries(${PROJECT_NAME} PUBLIC matrixio ${MPI_C_LIBRARIES} ${MKL_LIBRARIES_ILP})
