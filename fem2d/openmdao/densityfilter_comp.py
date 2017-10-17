import numpy as np

from openmdao.api import ExplicitComponent


class DensityFilterComp(ExplicitComponent):
    # this simplest density filter uses <8 elements (nearest neighbors)
    def initialize(self):
        self.metadata.declare('length_x', type_=(int, float), required=True)
        self.metadata.declare('length_y', type_=(int, float), required=True)
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)
        self.metadata.declare('num_dvs', type_=int, required=True)
        self.metadata.declare('radius', type_=float, required=True)
        
    def setup(self):
        num_dvs = self.metadata['num_dvs']
        length_x = self.metadata['length_x']
        length_y = self.metadata['length_y']
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']
        radius = self.metadata['radius']

        self.add_input('dvs', shape=num)
        self.add_output('dvs_bar', shape=num)

        num_elem = (num_nodes_x - 1) * (num_nodes_y - 1)

        for ii in range(0, num_nodes_x - 1):
            for jj in range(0, num_nodes_y - 1):
                


        
        

        self.declare_partials('dvs_bar', 'dvs', row= , col=)

    def compute(self, inputs, outputs):
        outputs['y'] = self.metadata['mtx'].dot(inputs['x'])
