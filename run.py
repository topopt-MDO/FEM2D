from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import scipy.sparse
import scipy.sparse.linalg as SL

from fem2d import PyFEMSolver

length_x = 6
length_y = 6
num_nodes_x = 6
num_nodes_y = 6
E = 1
nu = 0.5

a = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)

size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64 + 2 * 2 * num_nodes_y
data = np.zeros(size)
rows = np.zeros(size, np.int32)
cols = np.zeros(size, np.int32)

a.get_stiffness_matrix(data, rows, cols)

print(data)
print(rows)
print(cols)

orig_nodes = np.zeros((num_nodes_x, num_nodes_y, 2))
for ix in range(num_nodes_x):
    for iy in range(num_nodes_y):
        orig_nodes[ix, iy, 0] = ix / (num_nodes_x - 1) * length_x
        orig_nodes[ix, iy, 1] = iy / (num_nodes_y - 1) * length_y

def plot(axes):
    axes.clear()
    for ix in range(num_nodes_x - 1):
        for iy in range(num_nodes_y - 1):
            axes.plot(
                orig_nodes[ix+0:ix+2, iy+0, 0],
                orig_nodes[ix+0:ix+2, iy+0, 1],
                'k')
            axes.plot(
                orig_nodes[ix+0:ix+2, iy+1, 0],
                orig_nodes[ix+0:ix+2, iy+1, 1],
                'k')
            axes.plot(
                orig_nodes[ix+0, iy+0:iy+2, 0],
                orig_nodes[ix+0, iy+0:iy+2, 1],
                'k')
            axes.plot(
                orig_nodes[ix+1, iy+0:iy+2, 0],
                orig_nodes[ix+1, iy+0:iy+2, 1],
                'k')

size = num_nodes_x * num_nodes_y * 2 + 2 * num_nodes_y
mtx = scipy.sparse.csc_matrix((data, (rows, cols)), shape=(size, size))

# if 0:
axes = plt.gca()
plot(axes)
plt.show()

# if 1:
# plt.spy(mtx)
# plt.show()

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
    vecF = np.zeros(2*num_nodes_x*num_nodes_y + 2*num_nodes_y,);
    vecF[(num_nodes_y*(num_nodes_x-1)+int(num_nodes_y/2))*2+1] = - 0.01;
    return vecF

f = compute_force()

u = SL.spsolve(mtx,f)

print(u)
