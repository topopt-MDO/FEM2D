import numpy as np
import scipy.sparse 

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

        self.add_input('dvs', shape=num_dvs)
        self.add_output('dvs_bar', shape=num_dvs)

        num_elem_x = num_nodes_x - 1
        num_elem_y = num_nodes_y - 1
        num_elem = num_elem_x * num_elem_y

        lx = length_x / float(num_elem_x)
        ly = length_y / float(num_elem_y)

        cnt = 0

        row = []
        col = []
        wij = []
        
        cnt = 0
        for ix in range(0, num_elem_x):
            for iy in range(0, num_elem_y):
                eid_i  = ix * num_elem_y + iy
                x_i = (ix + 0.5) * lx
                y_i = (iy + 0.5) * ly
                dsum = 0
                dij_list = []
                negh_list = []
                num_neighbor = 0

                for jx in range(max(0, ix - 3), min(num_elem_x, ix + 3)):
                    for jy in range(max(0, iy - 3), min(num_elem_y, iy + 3)):
                        eid_j  = jx * num_elem_y + jy
                        x_j = (jx + 0.5) * lx
                        y_j = (jy + 0.5) * ly

                        dij = ((x_i-x_j)**2 + (y_i-y_j)**2)**0.5
                        # print(eid_i, eid_j, x_i, y_i, x_j, y_j, dij)
                        if dij < radius:
                            num_neighbor += 1
                            dij_list.append(radius - dij)
                            negh_list.append(eid_j)
                            dsum += radius - dij
                
                for tt in range(0, num_neighbor):
                    row.append(eid_i)
                    col.append(negh_list[tt])
                    wij.append(dij_list[tt]/dsum)
                    cnt += 1

        row = np.array(row)
        col = np.array(col)
        wij = np.array(wij)
        self.mtx = scipy.sparse.csr_matrix((wij, (row, col)), shape=(num_elem, num_elem))
                
        self.declare_partials('dvs_bar', 'dvs', rows = row, cols = col, val = wij)

    def compute(self, inputs, outputs):
        outputs['dvs_bar'] = self.mtx.dot(inputs['dvs'])
