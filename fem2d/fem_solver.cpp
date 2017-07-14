#include "fem_solver.h"
#include "LinAlg_MDO.cpp"
#include <assert.h>

FEMSolver::FEMSolver(
  int num_nodes_x, int num_nodes_y, double length_x, double length_y,
  double E, double nu
) {
  this->num_nodes_x = num_nodes_x;
  this->num_nodes_y = num_nodes_y;
  this->length_x = length_x;
  this->length_y = length_y;

  dr_dx = 2.0 / length_x * (num_nodes_x - 1);
  ds_dy = 2.0 / length_y * (num_nodes_y - 1);

  D_voigt.resize(3, vector<double>(3,0.0));

  compute_nodes();
  compute_elems();
  compute_D(E, nu);

  compute_Ke(1,2,3,4);
}

void FEMSolver::compute_D(double E, double nu){
  double sc = E/(1-nu*nu);

  D_voigt[0][0] = sc * 1;
  D_voigt[1][1] = sc * 1;
  D_voigt[0][1] = sc * nu;
  D_voigt[1][0] = sc * nu;
  D_voigt[2][2] = sc * 0.5 * (1 - nu);
}

FEMSolver::~FEMSolver() {
}

void FEMSolver::compute_elems() {
  int node_i, node_j;

  elems.resize(num_nodes_x, vector<vector<int> >(num_nodes_y, vector<int>(8, 0)));

  for (int i = 0; i < num_nodes_x - 1; i++){
    for (int j = 0; j < num_nodes_y - 1; j++){
      // -1, -1
      elems[i][j][0] = (i + 0) * num_nodes_y * 2 + (j + 0) * 2 + 0;
      elems[i][j][1] = (i + 0) * num_nodes_y * 2 + (j + 0) * 2 + 1;
      //  1, -1
      elems[i][j][2] = (i + 1) * num_nodes_y * 2 + (j + 0) * 2 + 0;
      elems[i][j][3] = (i + 1) * num_nodes_y * 2 + (j + 0) * 2 + 1;
      //  1,  1
      elems[i][j][4] = (i + 1) * num_nodes_y * 2 + (j + 1) * 2 + 0;
      elems[i][j][5] = (i + 1) * num_nodes_y * 2 + (j + 1) * 2 + 1;
      // -1,  1
      elems[i][j][6] = (i + 0) * num_nodes_y * 2 + (j + 1) * 2 + 0;
      elems[i][j][7] = (i + 0) * num_nodes_y * 2 + (j + 1) * 2 + 1;
    }
  }
}

void FEMSolver::compute_nodes() {
  nodes.resize(num_nodes_x, vector<vector<double> >(num_nodes_y, vector<double>(2,0))); // (x,y) for indx. (i,j)
  double sp_x = length_x/(num_nodes_x-1), sp_y = length_y/(num_nodes_y-1); // spacing btw. nodes
  for (int ii = 0; ii < num_nodes_x; ii++){
    for (int jj = 0; jj < num_nodes_y; jj++){
      nodes[ii][jj][0] = ii*sp_x;
      nodes[ii][jj][1] = jj*sp_y;
    }
  }

  // print(nodes);
}

void FEMSolver::compute_Ke(double x1, double x2, double y1, double y2){
  Matrix B;
  // N.resize(2, vector<double>(8,0.0));
  B.resize(3, vector<double>(8,0.0));

  // du_dx
  B[0][0] += -0.25 * dr_dx;
  B[0][1] += 0.;
  B[0][2] +=  0.25 * dr_dx;
  B[0][3] += 0.;
  B[0][4] +=  0.25 * dr_dx;
  B[0][5] += 0.;
  B[0][6] += -0.25 * dr_dx;
  B[0][7] += 0.;

  // dv_dy
  B[1][0] += 0;
  B[1][1] += -0.25 * ds_dy;
  B[1][2] += 0.;
  B[1][3] += -0.25 * ds_dy;
  B[1][4] += 0.;
  B[1][5] +=  0.25 * ds_dy;
  B[1][6] += 0.;
  B[1][7] +=  0.25 * ds_dy;

  // du_dy
  B[2][0] += -0.25 * ds_dy;
  B[2][1] += 0.;
  B[2][2] += -0.25 * ds_dy;
  B[2][3] += 0.;
  B[2][4] +=  0.25 * ds_dy;
  B[2][5] += 0.;
  B[2][6] +=  0.25 * ds_dy;
  B[2][7] += 0.;

  // dv_dx
  B[2][0] += 0;
  B[2][1] += -0.25 * dr_dx;
  B[2][2] += 0.;
  B[2][3] +=  0.25 * dr_dx;
  B[2][4] += 0.;
  B[2][5] +=  0.25 * dr_dx;
  B[2][6] += 0.;
  B[2][7] += -0.25 * dr_dx;

  Matrix BT = transpose(B);
  Matrix BD_tmp = dot(BT, D_voigt);
  Ke_ = dot(BD_tmp, B);
}

void FEMSolver::get_stiffness_matrix(double* data, int* rows, int* cols) {
  int index = 0;

  for (int ielem_x = 0; ielem_x < num_nodes_x - 1; ielem_x++) {
    for (int ielem_y = 0; ielem_y < num_nodes_y - 1; ielem_y++) {
      for (int imat_x = 0; imat_x < 8; imat_x++) {
        for (int imat_y = 0; imat_y < 8; imat_y++) {
          data[index] = Ke_[imat_x][imat_y];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = elems[ielem_x][ielem_y][imat_y];
          index += 1;
        }
      }
    }
  }

  // Lagrange multipliers
  int inode_x = 0;
  int idof = 0;
  int num_dofs = num_nodes_x * num_nodes_y * 2;

  for (int k = 0; k < 2; k++) {
    for (int inode_y = 0; inode_y < num_nodes_y; inode_y++) {
      idof = inode_x * num_nodes_y * 2 + inode_y * 2 + k;

      data[index] = 1.0;
      rows[index] = idof;
      cols[index] = idof + num_dofs;
      index += 1;

      data[index] = 1.0;
      rows[index] = idof + num_dofs;
      cols[index] = idof;
      index += 1;
    }
  }
  cout << index << "\n\n\n";
}

int main(){
  FEMSolver fem2(1,1,1,1,1,0.5);
  // for (int ii = 0; ii < 3; ii++){
  //   cout << "\n";
  //   for (int jj = 0; jj < 3; jj++){
  //     cout << fem2.D_voigt[ii][jj] << "\t";
  //   }
  //   cout << "\n";
  // }

}
