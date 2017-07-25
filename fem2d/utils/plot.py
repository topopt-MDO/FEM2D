import numpy as np
import matplotlib.pylab as plt


def _plot(nodes, axes, color='k'):
    for ix in range(nodes.shape[0] - 1):
        for iy in range(nodes.shape[1] - 1):
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


def get_mesh(num_nodes_x, num_nodes_y, length_x, length_y):
    orig_nodes = np.zeros((num_nodes_x, num_nodes_y, 2))
    for ix in range(num_nodes_x):
        for iy in range(num_nodes_y):
            orig_nodes[ix, iy, 0] = ix / (num_nodes_x - 1) * length_x
            orig_nodes[ix, iy, 1] = iy / (num_nodes_y - 1) * length_y

    return orig_nodes


def get_def_mesh(orig_nodes, disp, num_nodes_x, num_nodes_y, scale=1e0):
    return orig_nodes + disp.reshape((num_nodes_x, num_nodes_y, 2)) * scale


def plot_solution(orig_nodes, deflected_nodes=None):
    axes = plt.gca()
    _plot(orig_nodes, axes)
    if deflected_nodes is not None:
        _plot(deflected_nodes, axes, 'r')
    plt.show()
