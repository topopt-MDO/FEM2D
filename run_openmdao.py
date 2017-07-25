import numpy as np

from openmdao.api import Problem

from fem2d.fem2d import PyFEMSolver
from fem2d.openmdao.fem2d_group import FEM2DGroup
from fem2d.utils.plot import get_mesh, get_def_mesh, plot_solution
from fem2d.utils.forces import get_forces


num_nodes_x = 11
num_nodes_y = 11

length_x = 1
length_y = 1

E = 1000.
nu = 0.3
f = -10

fem_solver = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)
forces = get_forces(num_nodes_x, num_nodes_y, f=f)

model = FEM2DGroup(
    fem_solver=fem_solver, num_nodes_x=num_nodes_x, num_nodes_y=num_nodes_y, forces=forces)

prob = Problem(model)
prob.setup()
prob.run_model()

disp = prob['disp_comp.disp']

scale = 1e0
orig_nodes = get_mesh(num_nodes_x, num_nodes_y, length_x, length_y)
deflected_nodes = get_def_mesh(orig_nodes, disp, num_nodes_x, num_nodes_y, scale=scale)
plot_solution(orig_nodes, deflected_nodes=deflected_nodes)
