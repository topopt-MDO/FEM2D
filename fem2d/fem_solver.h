#include <iostream>
#include <vector>

using namespace std;

class FEMSolver {
public:
  FEMSolver(
    int num_nodes_x, int num_nodes_y, double length_x, double length_y,
    double E, double nu
  );
  ~FEMSolver();
  void get_stiffness_matrix(double* data, int* rows, int* cols);
private:
  int num_nodes, num_elems;
  int num_nodes_x, num_nodes_y;
  double length_x, length_y;
  vector<vector<double> > nodes;
  vector<vector<int> > elems;
  vector<vector<double> > Cijkl;
  void compute_nodes();
  void compute_elems();
  void compute_Cijkl(double E, double nu);
};
