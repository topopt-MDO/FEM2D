import numpy as np

from fem2d import PyFEMSolver

num_nodes_x = 1
num_nodes_y = 1 

a = PyFEMSolver(num_nodes_x, num_nodes_y,1,1,1,1)

size = num_nodes_x * num_nodes_y
#data = np.zeros(size)
#rows = np.zeros(size, int)
#cols = np.zeros(size, int)


data = np.zeros(size)
rows = np.zeros(size, dtype = 'int64')
cols = rows#np.ndarray(size, int)

print(data)
print(rows)
print(cols)


a.get_stiffness_matrix(data, rows, cols)

print(data)
print(rows)
print(cols)


