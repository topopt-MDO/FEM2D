import numpy as np

from fem2d import PyFEMSolver

num_nodes_x = 2
num_nodes_y = 3

a = PyFEMSolver(num_nodes_x, num_nodes_y,1,1,1,1)

size = num_nodes_x * num_nodes_y
data = np.zeros(size)
rows = np.zeros(size, np.int32)
cols = np.zeros(size, np.int32)

a.get_stiffness_matrix(data, rows, cols)

print(data)
print(rows)
print(cols)
