import numpy as np

from fem2d import PyFEMSolver

num_nodes_x = 3
num_nodes_y = 3

a = PyFEMSolver(num_nodes_x, num_nodes_y,1,1,1,0.5)

size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64 + 2 * 2 * num_nodes_y
data = np.zeros(size)
rows = np.zeros(size, np.int32)
cols = np.zeros(size, np.int32)

a.get_stiffness_matrix(data, rows, cols)

print(data)
print(rows)
print(cols)

import scipy.sparse
import matplotlib.pylab as plt
size = num_nodes_x * num_nodes_y * 2 + 2 * num_nodes_y
mtx = scipy.sparse.csc_matrix((data, (rows, cols)), shape=(size, size))
plt.spy(mtx)
plt.show()

# print(data)
# print(rows)
# print(cols)
