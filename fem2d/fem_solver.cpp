#include "fem_solver.h"

FEMSolver::FEMSolver(
  int num_nodes_x, int num_nodes_y, double length_x, double length_y,
  double E, double nu
) {
  this->num_nodes_x = num_nodes_x;
  this->num_nodes_y = num_nodes_y;
  this->length_x = length_x;
  this->length_y = length_y;

  compute_nodes();
  compute_elems();
  compute_Cijkl(E, nu);
}

FEMSolver::~FEMSolver() {
}

void FEMSolver::compute_elems() {}

void FEMSolver::compute_nodes() {}

void FEMSolver::compute_Cijkl(double E, double nu) {}

void FEMSolver::get_stiffness_matrix(double* data, int* rows, int* cols) {
  int size = num_nodes_x * num_nodes_y;

  for (int i = 0; i < size; i++){
    data[i] = 1;
    rows[i] = 2;
    cols[i] = 3;
  }
}
