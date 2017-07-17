from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import scipy.sparse
import scipy.sparse.linalg
import scipy.linalg

from fem2d import PyFEMSolver

length_x = 1
length_y = 1
num_nodes_x = 10
num_nodes_y = 10
E = 1
nu = 0

a = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)

size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64 + 2 * 2 * num_nodes_y
data = np.zeros(size)
rows = np.zeros(size, np.int32)
cols = np.zeros(size, np.int32)

a.get_stiffness_matrix(data, rows, cols)

orig_nodes = np.zeros((num_nodes_x, num_nodes_y, 2))
for ix in range(num_nodes_x):
    for iy in range(num_nodes_y):
        orig_nodes[ix, iy, 0] = ix / (num_nodes_x - 1) * length_x
        orig_nodes[ix, iy, 1] = iy / (num_nodes_y - 1) * length_y

def plot(nodes, axes, color='k'):
    for ix in range(num_nodes_x - 1):
        for iy in range(num_nodes_y - 1):
            axes.plot(
                nodes[ix+0:ix+2, iy+0, 0],
                nodes[ix+0:ix+2, iy+0, 1],
                color)
            axes.plot(
                nodes[ix+0:ix+2, iy+1, 0],
                nodes[ix+0:ix+2, iy+1, 1],
                color)
            axes.plot(
                nodes[ix+0, iy+0:iy+2, 0],
                nodes[ix+0, iy+0:iy+2, 1],
                color)
            axes.plot(
                nodes[ix+1, iy+0:iy+2, 0],
                nodes[ix+1, iy+0:iy+2, 1],
                color)

size = num_nodes_x * num_nodes_y * 2 + 2 * num_nodes_y
mtx = scipy.sparse.csc_matrix((data, (rows, cols)), shape=(size, size))

<<<<<<< HEAD
# if 0:
axes = plt.gca()
plot(axes)
plt.show()

# if 1:
# plt.spy(mtx)
# plt.show()
=======
# print(np.min(np.unique(rows)), np.max(np.unique(rows)))
# print(np.min(np.unique(cols)), np.max(np.unique(cols)))
print(size)

print('x',np.linalg.matrix_rank(mtx.todense()[:8,:8]))

# exit()

if 0:
    axes = plt.gca()
    plot(orig_nodes, axes)
    plt.show()

if 0:
    plt.spy(mtx)
    plt.show()
>>>>>>> 2a1e9c36480d580c259ee2a0a3c8af2b9a9ae558

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
    vecF[(num_nodes_y*(num_nodes_x-1)+int(num_nodes_y/2))*2+1] = -1e0;
    return vecF

f = compute_force()

# print(mtx.todense()[:8,:8])
# print(mtx.todense())

if 1:
    lu = scipy.sparse.linalg.splu(mtx)
    u = lu.solve(f)
else:
    lu = scipy.linalg.lu_factor(mtx.todense())
    u = scipy.linalg.lu_solve(lu, f)

deflected_nodes = u[:num_nodes_x*num_nodes_y*2].reshape((num_nodes_x, num_nodes_y, 2))

if 1:
    axes = plt.gca()
    plot(orig_nodes, axes)
    plot(orig_nodes + deflected_nodes * 1e0, axes, 'r')
    plt.show()
