import numpy as np

from openmdao.api import ExplicitComponent

from fem2d.utils.plot import get_mesh, plot_imshow


class PenalizationComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('p', type_=(int, float), required=True)
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)
        self.metadata.declare('nodes', type_=np.ndarray, required=True)

    def setup(self):
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        self.mesh = self.metadata['nodes']
        self.counter = 0

        multiplier_size = (num_nodes_x - 1) * (num_nodes_y - 1)

        self.add_input('densities', shape=multiplier_size)
        self.add_output('multipliers', shape=multiplier_size)

        arange = np.arange(multiplier_size)
        self.declare_partials('multipliers', 'densities', rows=arange, cols=arange)

    def compute(self, inputs, outputs):
        p = self.metadata['p']

        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']
        if self.counter % 5 == 0:
            plot_imshow(
                self.mesh, inputs['densities'].reshape((num_nodes_x - 1, num_nodes_y - 1)),
                save='save/save%03i.png'%self.counter)
        self.counter += 1

        outputs['multipliers'] = inputs['densities'] ** p

    def compute_partials(self, inputs, outputs, partials):
        p = self.metadata['p']

        partials['multipliers', 'densities'] = p ** inputs['densities'] ** (p - 1)
