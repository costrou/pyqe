"""Cell Namelist

"""
from pyqe.namelist import Namelist

# Functions to evaluate range
def isPositive(value):
    return value > 0.0

class Cell(Namelist):
    """Cell Namelist

    """
    def _defaultCellDynamics(self):
        pass #TODO needs qe.control to set default

    def __init__(self):
        name = "CELL"
        keypairs = {}
        keys = {
            'cell_dynamics': [0, str, self._defaultCellDynamics, ('none', 'sd', 'damp-pr', 'bfgs', 'pr', 'w'), None],
            'press': [0, float, 0.0, None, None],
            'wmass': [0, float, None, isPositive, None], #TODO Check onlin if vc-md or vd-relax
            'cell_factor': [0, float, 1.2, None, None],
            'press_conv_thr': [0, float, 0.5, None, None],
            'cell_dofree': [0, str, 'all', ('all', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape'), None]
        }
        super().__init__(name, keypairs, keys)
