import numpy as np

from openmdao.api import Group, IndepVarComp

from fem2d.openmdao.states_comp import StatesComp
from fem2d.openmdao.disp_comp import DispComp
from fem2d.openmdao.compliance_comp import ComplianceComp
from fem2d.fem2d import PyFEMSolver


class FEM2DGroup(Group):

    def initialize(self):
        self.metadata.declare('fem_solver', type_=PyFEMSolver, required=True)
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)
        self.metadata.declare('forces', type_=np.ndarray, required=True)

    def setup(self):
        fem_solver = self.metadata['fem_solver']
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']
        forces = self.metadata['forces']

        state_size = 2 * num_nodes_x * num_nodes_y + 2 * num_nodes_y
        disp_size = 2 * num_nodes_x * num_nodes_y

        rhs = np.zeros(state_size)
        rhs[:disp_size] = forces

        # inputs
        comp = IndepVarComp()
        comp.add_output('rhs', val=rhs)
        comp.add_output('forces', val=forces)
        self.add_subsystem('inputs_comp', comp)

        # states
        comp = StatesComp(fem_solver=fem_solver, num_nodes_x=num_nodes_x, num_nodes_y=num_nodes_y)
        self.add_subsystem('states_comp', comp)
        self.connect('inputs_comp.rhs', 'states_comp.rhs')

        # disp
        comp = DispComp(num_nodes_x=num_nodes_x, num_nodes_y=num_nodes_y)
        self.add_subsystem('disp_comp', comp)
        self.connect('states_comp.states', 'disp_comp.states')

        # compliance
        comp = ComplianceComp(num_nodes_x=num_nodes_x, num_nodes_y=num_nodes_y)
        self.add_subsystem('compliance_comp', comp)
        self.connect('disp_comp.disp', 'compliance_comp.disp')
        self.connect('inputs_comp.forces', 'compliance_comp.forces')
