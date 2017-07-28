import numpy as np

from openmdao.api import ExplicitComponent

from fem2d.utils.gauss_quadrature import gauss_wts_dict


class AveragingComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)
        self.metadata.declare('quad_order', type_=int, required=True)

    def setup(self):
        nx = self.metadata['num_nodes_x'] - 1
        ny = self.metadata['num_nodes_y'] - 1
        quad_order = self.metadata['quad_order']

        self.gauss_wts = gauss_wts = gauss_wts_dict[quad_order] / 2.

        self.add_input('x', shape=nx * quad_order * ny * quad_order)
        self.add_output('y', shape=nx * ny)

        data = np.einsum('ik,j,l->ijkl',
            np.ones((nx, ny)), gauss_wts, gauss_wts,
        ).flatten()
        rows = np.einsum('ik,jl->ijkl',
            np.arange(nx * ny).reshape((nx, ny)), np.ones((quad_order, quad_order)),
        ).flatten()
        cols = np.arange(nx * quad_order * ny * quad_order)
        self.declare_partials('y', 'x', val=data, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        nx = self.metadata['num_nodes_x'] - 1
        ny = self.metadata['num_nodes_y'] - 1
        quad_order = self.metadata['quad_order']

        outputs['y'] = np.einsum('ijkl,j,l->ik',
            inputs['x'].reshape((nx, quad_order, ny, quad_order)), self.gauss_wts, self.gauss_wts,
        ).flatten()
