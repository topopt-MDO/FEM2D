from libcpp.vector cimport vector
import numpy as np
cimport numpy as np


cdef extern from "fem_solver.h":
  cdef cppclass FEMSolver:
    FEMSolver(int, int, double, double, double, double) except +
    void get_stiffness_matrix(double* data, int* rows, int* cols)

cdef class PyFEMSolver:

    cdef FEMSolver *thisptr
    def __cinit__(self, int num_nodes_x, int num_nodes_y, double length_x, double length_y,
            double E, double nu,
        ):
        self.thisptr = new FEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)
    def __dealloc__(self):
        del self.thisptr
    def get_stiffness_matrix(
            self, np.ndarray[double] data, np.ndarray[int] rows, np.ndarray[int] cols):
        self.thisptr.get_stiffness_matrix(&data[0], &rows[0], &cols[0])
