"""Electrons Namelist

"""
from pyqe.namelist import Namelist


# Functions to evaluate range
def isPositive(value):
    return value > 0.0

class Electrons(Namelist):
    """Electrons Namelist

    """
    def _defaultStartingpot(self):
        pass #TODO needs control.calculation 
    
    def __init__(self):
        name = "ELECTRONS"
        keypairs = {}
        keys = {
            'electron_maxstep': [0, int, 100, isPositive, None],
            'scf_must_converge': [0, bool, True, None, None],
            'conv_thr': [0, float, 1e-6, isPositive, None],
            'adaptive_thr': [0, bool, False, None, None],
            'conv_thr_init': [0, float, 1e-3, isPositive, None],
            'conv_thr_multi': [0, float, 1e-1, isPositive, None],
            'mixing_mode': [0, str, 'plain', ('plain', 'TF', 'local-TF'), None],
            'mixing_beta': [0, float, 0.7, None, None],
            'mixing_ndim': [0, int, 8, isPositive, None],
            'mixing_fixed_ns': [0, int, 0, isPositive, None],
            'diagonalization': [0 , str, 'david', ('david', 'cg', 'cg-serial'), None],
            'ortho_para': [0, int, 0, None, None],
            'diago_thr_init': [0, float, None, None, None],
            'diago_cg_maxiter': [0, int, None, None, None],
            'diago_david_ndim': [0, int, 4, None, None],
            'diago_full_acc': [0, bool, False, None, None],
            'efield': [0, float, 0.0, None, None],
            'efield_cart': [1, float, 0.0, None, None],
            'startingpot': [0, str, self._defaultStartingpot, ('atomic', 'file'), None],
            'startingwfc': [0, str, 'atomic+random', ('atomic', 'atomic+random', 'random', 'file'), None],
            'tqr': [0, bool, False, None, None],
        }
        Namelist.__init__(self, name, keys)
