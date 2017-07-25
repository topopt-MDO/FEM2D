import numpy as np
import scipy.sparse
import scipy.sparse.linalg

from openmdao.api import ImplicitComponent

from fem2d.fem2d import PyFEMSolver


class StatesComp(ImplicitComponent):

    def initialize(self):
        self.metadata.declare('fem_solver', type_=PyFEMSolver, required=True)
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)

    def setup(self):
        fem_solver = self.metadata['fem_solver']
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        state_size = 2 * num_nodes_x * num_nodes_y + 2 * num_nodes_y

        self.add_input('rhs', shape=state_size)
        self.add_output('states', shape=state_size)

        size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64 + 2 * 2 * num_nodes_y
        self.data = data = np.zeros(size)
        self.rows = rows = np.zeros(size, np.int32)
        self.cols = cols = np.zeros(size, np.int32)

        fem_solver.get_stiffness_matrix(data, rows, cols)
        self.declare_partials('states', 'states', rows=rows, cols=cols)

        data = -np.ones(state_size)
        arange = np.arange(state_size)
        self.declare_partials('states', 'rhs', val=data, rows=arange, cols=arange)

    def _get_mtx(self, inputs):
        fem_solver = self.metadata['fem_solver']
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        state_size = 2 * num_nodes_x * num_nodes_y + 2 * num_nodes_y

        data, rows, cols = self.data, self.rows, self.cols

        fem_solver.get_stiffness_matrix(data, rows, cols)

        mtx = scipy.sparse.csc_matrix((data, (rows, cols)), shape=(state_size, state_size))

        return mtx

    def apply_nonlinear(self, inputs, outputs, residuals):
        self.mtx = self._get_mtx(inputs)

        residuals['states'] = self.mtx.dot(outputs['states']) - inputs['rhs']

    def solve_nonlinear(self, inputs, outputs):
        self.mtx = self._get_mtx(inputs)
        self.lu = scipy.sparse.linalg.splu(self.mtx)

        outputs['states'] = self.lu.solve(inputs['rhs'])

    def linearize(self, inputs, outputs, partials):
        partials['states', 'states'] = self.data

    def solve_linear(self, d_outputs, d_residuals, mode):
        if mode == 'fwd':
            d_outputs['states'] = self.lu.solve(d_residuals['states'], 'N')
        elif mode == 'rev':
            d_residuals['states'] = self.lu.solve(d_outputs['states'], 'T')
