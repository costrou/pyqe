"""Ions Namelist

"""
from pyqe.namelist import Namelist

# Functions to evaluate range
def isPositive(value):
    return value > 0.0


class Ions(Namelist):
    """Ions Namelist

    """
    def _defaultIonDynamics():
        pass #TODO

    def _checkPotExtrapolation(self, qe):
        return [True, None] #TODO

    def _checkWfcExtrapolation(self, qe):
        return [True, None] #TODO
    
    def __init__(self):
        name = "IONS"
        keys = {
            'ion_dynamics': [0, str, self._defaultIonDynamics, ('bfgs', 'damp', 'verlet', 'langevin', 'langevin-smc', 'beeman'), None],
            'ion_positions': [0, str, 'default', ('default', 'from_input'), None],
            'pot_extrapolation': [0, str, 'atomic', ('none', 'atomic', 'first_order', 'second_order'), self._checkPotExtrapolation], 
            'wfc_extrapolation': [0, str, 'none', ('none', 'first-order', 'second-order'), self._checkWfcExtrapolation],
            'remove_rigid_rot': [0, bool, False, (), None],
            'ion_temperature': [0, str, 'not-controlled', ('rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'berendsen', 'andersen', 'initial', 'not_controlled'), None],
            'tempw': [0, float, 300.0, isPositive, None],
            'tolp': [0, float, 100.0, isPositive, None],
            'delta_t': [0, float, 1.0, None, None],
            'nraise': [0, int, 1, None, None],
            'refold_pos': [0, bool, False, None, None],
            'upscale': [0, float, 100.0, None, None], #TODO Check only in bfgs calc
            'bfgs_ndim': [0, int, 1, isPositive, None], #TODO Check only in bfgs calc
            'trust_radius_max': [0, float, 0.8, isPositive, None], #TODO Check only in bfgs calc
            'trust_radius_min': [0, float, 1e-3, isPositive, None], #TODO Check only in bfgs calc
            'trust_radius_ini': [0, float, 0.5, None, None], #TODO Check only in bfgs calc
            'w_1': [0, float, 0.01, None, None], #TODO Check only in bfgs calc
            'w_2': [0, float, 0.5, None, None] #TODO Check only in bfgs calc
        }
        Namelist.__init__(self, name, keys)
