#pragma once
#ifndef SLEPCTEST_SOLVERS_TESTS
#define SLEPCTEST_SOLVERS_TESTS

#include <mkl.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <cassert>

using namespace std;
namespace cr = std::chrono;

enum class MatrixFileFormat {
  CSR,
  MTX
};

class SolverTest {
public:
  int argc;
  char** argv;

  string A_name;
  MatrixFileFormat file_format;

  int maxiter = 1000;
  double tol = 1e-8;

  virtual int test() = 0;

protected:
  MKL_INT sizea = 0;
  MKL_INT nnza = 0;
  vector<MKL_INT> ia, ja;
  vector<double> a;

  cr::time_point<cr::system_clock> start, end;

  SolverTest() {};

  void read_matrix_file(const std::string& file_name,
    const MatrixFileFormat file_format,
    std::vector<MKL_INT>& ia, std::vector<MKL_INT>& ja,
    std::vector<double>& A, MKL_INT& matrix_size, MKL_INT& nnz);
};

class LAESTest : public SolverTest {
public:
  virtual int test();

protected:
  LAESTest() {};
  virtual int testSpecific() = 0;
};

class PETScCGTest : public LAESTest {
protected:
  string petscfilename;
  cr::time_point<cr::system_clock> scatter_start, scatter_end;

  int test();
  int testSpecific();
};

#endif