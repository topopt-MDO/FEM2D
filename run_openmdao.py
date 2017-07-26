import numpy as np

from openmdao.api import Problem, view_model, ScipyOptimizer, pyOptSparseDriver

from fem2d.fem2d import PyFEMSolver
from fem2d.openmdao.fem2d_group import FEM2DGroup
from fem2d.utils.plot import get_mesh, plot_solution, plot_contour, plot_imshow
from fem2d.utils.forces import get_forces


num_nodes_x = 41
num_nodes_y = 21

length_x = 2
length_y = 1

E = 1000.
nu = 0.3
f = -10
p = 3
w = 0.9

fem_solver = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)

forces = get_forces(num_nodes_x, num_nodes_y, f=f)
nodes = get_mesh(num_nodes_x, num_nodes_y, length_x, length_y)

model = FEM2DGroup(
    fem_solver=fem_solver, num_nodes_x=num_nodes_x, num_nodes_y=num_nodes_y, forces=forces, p=p,
    nodes=nodes, w=w)

prob = Problem(model)

# prob.driver = ScipyOptimizer()
# prob.driver.options['optimizer'] = 'SLSQP'
# prob.driver.options['tol'] = 1e-9
# prob.driver.options['disp'] = True

prob.driver = pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings['Major optimality tolerance'] = 1e-6
prob.driver.opt_settings['Major feasibility tolerance'] = 1e-7
prob.driver.opt_settings['Verify level'] = -1

prob.setup()

prob.run_driver()
# prob.run_model()
# prob.check_partials(compact_print=True)

# view_model(prob)

disp = prob['disp_comp.disp']
densities = prob['inputs_comp.densities'].reshape((num_nodes_x - 1, num_nodes_y - 1))

nodal_densities = np.zeros((num_nodes_x, num_nodes_y))
nodal_densities[:-1, :-1] += densities
nodal_densities[ 1:, :-1] += densities
nodal_densities[:-1,  1:] += densities
nodal_densities[ 1:,  1:] += densities
nodal_densities[1:-1, 1:-1] /= 4.
nodal_densities[1:-1,  0] /= 2.
nodal_densities[1:-1, -1] /= 2.
nodal_densities[ 0, 1:-1] /= 2.
nodal_densities[-1, 1:-1] /= 2.

scale = 1e0
disp2 = disp.reshape((num_nodes_x, num_nodes_y, 2))[:, :, 1] * scale
# plot_solution(orig_nodes, deflected_nodes=deflected_nodes)
# plot_contour(nodes, field=nodal_densities)
plot_imshow(nodes, densities)
