import numpy as np
import numpy.linalg as LA

from fem2d import PyFEMSolver

num_nodes_x = 4
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

# def compute_nodes(num_nodes_x, num_nodes_y, length_x, length_y):
#     nodes = np.zeros([num_nodes_x, num_nodes_y, 2])
#     sp_x = length_x/(num_nodes_x-1)
#     sp_y = length_y/(num_nodes_y-1)
#     idx = np.linspace(0,length_x,num_nodes_x)
#     idy = np.linspace(0,length_y,num_nodes_y)
#     [x,y] = np.meshgrid(idx,idy)

#     print(x)
#     print(y)
    
#     pass

# compute_nodes(2,3,1,4)



def compute_force():
    vecF = np.zeros(2*num_nodes_x*num_nodes_y,);
    vecF[(num_nodes_y*(num_nodes_x-1)+int(num_nodes_y/2))*2+1] = - 5;
    return vecF

f = compute_force()

