import numpy as np


def get_rbf_mtx(num_out_x, num_out_y, num_in_x, num_in_y):
    num_out = num_out_x * num_out_y
    num_in = num_in_x * num_in_y

    # ----- ----- ----- ----- ----- ----- -----

    coord_nodes = np.zeros((num_out_x, num_out_y, 2))
    coord_param = np.zeros((num_in_x, num_in_y, 2))

    coord_nodes[:, :, 0] = np.einsum('i,j->ij',
        np.linspace(0, 1, num_out_x), np.ones(num_out_y))
    coord_nodes[:, :, 1] = np.einsum('i,j->ij',
        np.ones(num_out_x), np.linspace(0, 1, num_out_y))

    coord_param[:, :, 0] = np.einsum('i,j->ij',
        np.linspace(0, 1, num_in_x), np.ones(num_in_y))
    coord_param[:, :, 1] = np.einsum('i,j->ij',
        np.ones(num_in_x), np.linspace(0, 1, num_in_y))

    coord_nodes = coord_nodes.reshape((num_out_x * num_out_y, 2))
    coord_param = coord_param.reshape((num_in_x * num_in_y, 2))

    # ----- ----- ----- ----- ----- ----- -----

    grid_nodes = np.zeros((num_out, num_in, 2))
    grid_param = np.zeros((num_out, num_in, 2))

    for ind in range(2):
        grid_nodes[:, :, ind] = np.einsum('i,j->ij', coord_nodes[:, ind], np.ones(num_in))
        grid_param[:, :, ind] = np.einsum('i,j->ij', np.ones(num_out), coord_param[:, ind])

    diff_sq = (grid_nodes - grid_param) ** 2
    r_sq = diff_sq[:, :, 0] + diff_sq[:, :, 1]
    phi_mtx = np.exp(-r_sq)

    return phi_mtx
