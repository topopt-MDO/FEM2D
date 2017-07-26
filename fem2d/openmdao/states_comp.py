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
        multiplier_size = (num_nodes_x - 1) * (num_nodes_y - 1)

        self.add_input('multipliers', shape=multiplier_size)
        self.add_input('rhs', shape=state_size)
        self.add_output('states', shape=state_size)

        size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64 + 2 * 2 * num_nodes_y
        self.data = data = np.zeros(size)
        self.rows = rows = np.zeros(size, np.int32)
        self.cols = cols = np.zeros(size, np.int32)

        fem_solver.get_stiffness_matrix(np.ones(multiplier_size), data, rows, cols)
        self.declare_partials('states', 'states', rows=rows, cols=cols)

        size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64
        self.data_d = data_d = np.zeros(size)
        self.rows_d = rows_d = np.zeros(size, np.int32)
        self.cols_d = cols_d = np.zeros(size, np.int32)

        fem_solver.get_stiffness_matrix_derivs(np.ones(state_size), data_d, rows_d, cols_d)
        self.declare_partials('states', 'multipliers', rows=rows_d, cols=cols_d)

        data = -np.ones(state_size)
        arange = np.arange(state_size)
        self.declare_partials('states', 'rhs', val=data, rows=arange, cols=arange)

    def _get_mtx(self, inputs):
        fem_solver = self.metadata['fem_solver']
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        state_size = 2 * num_nodes_x * num_nodes_y + 2 * num_nodes_y

        data, rows, cols = self.data, self.rows, self.cols

        fem_solver.get_stiffness_matrix(inputs['multipliers'], data, rows, cols)

        mtx = scipy.sparse.csc_matrix((data, (rows, cols)), shape=(state_size, state_size))

        return mtx

    def _compute_mtx_derivs(self, outputs):
        fem_solver = self.metadata['fem_solver']

        data_d, rows_d, cols_d = self.data_d, self.rows_d, self.cols_d

        fem_solver.get_stiffness_matrix_derivs(outputs['states'], data_d, rows_d, cols_d)

    def _solve(self, sol, rhs, mode):
        if mode == 'fwd':
            arg = 'N'
        elif mode == 'rev':
            arg = 'T'

        class Callback(object):
            def __init__(self, mtx):
                self.counter = 0
                self.mtx = mtx
            def __call__(self, xk):
                # print('%3i ' % self.counter, np.linalg.norm(self.mtx.dot(xk) - rhs))
                print('%3i ' % self.counter, np.linalg.norm(xk))
                self.counter += 1

        class PC(object):
            def __init__(self, arg, ilu):
                self.arg = arg
                self.ilu = ilu
            def __call__(self, rhs):
                return self.ilu.solve(rhs, self.arg)

        size = sol.shape[0]
        pc_op = scipy.sparse.linalg.LinearOperator((size, size), matvec=PC(arg, self.ilu))
        sol[:] = scipy.sparse.linalg.gmres(
            self.mtx.T, rhs, x0=sol, M=pc_op,
            callback=Callback(self.mtx), tol=1e-10, restart=200,
        )[0]

        # sol[:] = self.ilu.solve(rhs, arg)

    def apply_nonlinear(self, inputs, outputs, residuals):
        self.mtx = self._get_mtx(inputs)

        residuals['states'] = self.mtx.dot(outputs['states']) - inputs['rhs']

    def solve_nonlinear(self, inputs, outputs):
        self.mtx = self._get_mtx(inputs)
        self.ilu = scipy.sparse.linalg.spilu(self.mtx, drop_tol=1e-10)

        self._solve(outputs['states'], inputs['rhs'], 'fwd')
        self._compute_mtx_derivs(outputs)

    def linearize(self, inputs, outputs, partials):
        self.mtx = self._get_mtx(inputs)
        self._compute_mtx_derivs(outputs)

        partials['states', 'states'] = self.data
        partials['states', 'multipliers'] = self.data_d

        self.ilu = scipy.sparse.linalg.spilu(self.mtx, drop_tol=1e-10)

    def solve_linear(self, d_outputs, d_residuals, mode):
        if mode == 'fwd':
            self._solve(d_outputs['states'], d_residuals['states'], 'fwd')
        elif mode == 'rev':
            self._solve(d_residuals['states'], d_outputs['states'], 'rev')
