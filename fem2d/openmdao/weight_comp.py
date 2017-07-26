import numpy as np

from openmdao.api import ExplicitComponent


class WeightComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)

    def setup(self):
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        multiplier_size = (num_nodes_x - 1) * (num_nodes_y - 1)

        self.add_input('multipliers', shape=multiplier_size)
        self.add_output('weight')
        self.declare_partials('weight', 'multipliers', val=np.ones((1, multiplier_size)))

    def compute(self, inputs, outputs):
        outputs['weight'] = sum(inputs['multipliers'])
