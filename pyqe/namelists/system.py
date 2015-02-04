"""System Namelist

"""
from pyqe.namelist import Namelist

# Functions to evaluate range
def isPositive(value):
    return value > 0.0

def isWithinOneOfZero(value):
    return abs(value) <= 1.0

def isBtwZeroOne(value):
    return value >= 0.0 and value <= 1.0

def isGTFour(value):
    return value > 4

class System(Namelist):
    """System Namelist

    """

    def _rangeIbrav(self, value):
        return value in ([i for i in range(15)] + [-5])

    def _rangeSpaceGroup(self, value):
        return value in (i for i in range(231))

    def _defaultEcutrho(self):
        return 4.0 * self.get_current_value('ecutwfc')

    def _defaultEcutfock(self):
        return self.get_current_value('ecutrho')

    def _defaultnqx1(self):
        return self.get_current_value('nr1')

    def _defaultnqx2(self):
        return self.get_current_value('nr2')

    def _defaultnqx3(self):
        return self.get_current_value('nr3')

    def _checkIbrav(self, qe):
        return [True, None]

    def _checkCelldm(self, qe):
        return [True, None]

    def _checkA(self, qe):
        return [True, None]

    def _checkB(self, qe):
        return [True, None]

    def _checkC(self, qe):
        return [True, None]

    def _checkCosAB(self, qe):
        return [True, None]

    def _checkCosAC(self, qe):
        return [True, None]

    def _checkCosBC(self, qe):
        return [True, None]

    def _checkNat(self, qe):
        return [True, None]

    def _checkNtyp(self, qe):
        return [True, None]

    def _checkNri(self, qe):
        if not qe.control.get_set_value('nr1') and \
           qe.control.get_set_value('nr2') and \
            qe.control.get_set_value('nr3'):
            return [False, "Must set 'nr1', 'nr2', and 'nr3'"]
        return [True, None]

    def _checkNris(self, qe):
        if not qe.control.get_set_value('nr1s') and \
           qe.control.get_set_value('nr2s') and \
            qe.control.get_set_value('nr3s'):
            return [False, "Must set 'nr1s', 'nr2s', and 'nr3s'"]
        return [True, None]

    def _checkSpaceGroup(self, qe):
        return [True, None]

    def __init__(self):
        name = "SYSTEM"
        keys = {
            'ibrav': [0, int, None, self._rangeIbrav, self._checkIbrav],
            'celldm': [1, float, None, isPositive, self._checkCelldm],
            'A': [0, float, None, isPositive, self._checkA],
            'B': [0, float, None, isPositive, self._checkB],
            'C': [0, float, None, isPositive, self._checkC],
            'cosAB': [0, float, None, isWithinOneOfZero, self._checkCosAB],
            'cosAC': [0, float, None, isWithinOneOfZero, self._checkCosAC],
            'cosBC': [0, float, None, isWithinOneOfZero, self._checkCosBC],
            'nat': [0, int, None, isPositive, self._checkNat],
            'ntyp': [0, int, None, isPositive, self._checkNtyp],
            'nbnd': [0, int, None, isGTFour, None], #TODO default nbnd (insulator)
            'tot_charge': [0, float, 0.0, None, None],
            'tot_magnetization': [0, float, -1.0, isWithinOneOfZero, None],
            'starting_magnetization': [1, float, None, isWithinOneOfZero, None],
            'ecutwfc': [0, float, None, isPositive, None],
            'ecutrho': [0, float, self._defaultEcutrho, isPositive, None],
            'ecutfock': [0, float, self._defaultEcutfock, isPositive, None],
            'nr1': [0, int, None, isPositive, self._checkNri],
            'nr2': [0, int, None, isPositive, self._checkNri],
            'nr3': [0, int, None, isPositive, self._checkNri],
            'nr1s': [0, int, None, isPositive, self._checkNris],
            'nr2s': [0, int, None, isPositive, self._checkNris],
            'nr3s': [0, int, None, isPositive, self._checkNris],
            'nosym': [0, bool, False, None, None],
            'nosym_evc': [0, bool, False, None, None],
            'noinv': [0, bool, False, None, None],
            'no_t_rev': [0, bool, False, None, None],
            'force_symmorphic': [0, bool, False, None, None],
            'use_all_frac': [0, bool, False, None, None],
            'occupations': [0, str, None, ('smearing', 'tetrahedra', 'fixed', 'from_input'), None],
            'one_atom_occupations': [0, bool, False, None, None],
            'starting_spin_angle': [0, bool, False, None, None],
            'degauss': [0, float, 0.0, isPositive, None],
            'smearing': [0, str, 'gaussian', ('gaussian', 'methfessel-paxton', 'm-p', 'mp', 'mazari-vanderbilt', 'cold', 'm-v', 'mv', 'fermi-dirac', 'f-d', 'fd'), None],
            'nspin': [0, int, 1, (1, 2, 4), None],
            'noncolin': [0, bool, False, None, None],
            'ecfixed': [0, float, 0.0, isPositive, None],
            'qcutz': [0, float, 0.0, isPositive, None],
            'q2sigma': [0, float, 0.1, isPositive, None],
            'input_dft': [0, str, None, None, None],
            'exx_fraction': [0, float, None, isBtwZeroOne, None],
            'screening_parameter': [0, float, 0.106, None, None],
            'exxdiv_treatment': [0, str, 'gygi-baldereshi', ('gygi-baldereschi', 'vcut_spherical', 'vcut_ws', 'none'), None],
            'x_gamma_extrapolation': [0, bool, True, None, None],
            'ecutvcut': [0, float, 0.0, isPositive, None],
            'nqx1': [0, int, self._defaultnqx1, isPositive, None],
            'nqx2': [0, int, self._defaultnqx2, isPositive, None],
            'nqx3': [0, int, self._defaultnqx3, isPositive, None],
            'lda_plus_u': [0, bool, False, None, None],
            'lda_plus_u_kind': [0, int, 0, (0, 1), None],
            'Hubbard_U': [1, float, 0.0, None, None],  #TODO valid range TODO check
            'Hubbard_J0': [1, float, 0.0, None, None], #TODO valid range TODO check
            'Hubbard_alpha': [1, float, 0.0, None, None], #TODO valid range TODO check
            'Hubbard_beta': [1, float, 0.0, None, None], #TODO valid range TODO check
            'Hubbard_J': [2, float, 0.0, None, None], #TODO valid range TODO check
            'starting_ns_eigenvalue': [3, float, -1.0, None, None], #TODO valid range TODO check
            'U_projection_type': [0, str, 'atomic', ('atomic', 'ortho-atomic', 'norm-atomic', 'file', 'pseduo'), None],
            'edir': [0, int, None, (1, 2, 3), None],
            'emaxpos': [0, float, 0.5, isBtwZeroOne, None],
            'eopreg': [0, float, 0.1, isBtwZeroOne, None], 
            'eamp': [0, float, 0.001, None, None], #TODO range
            'angle1': [1, float, None, None, None], #TODO check (1 .. ntyp)
            'angle2': [1, float, None, None, None], #TODO check (1 .. ntyp)
            'constrained_magnetization': [0, str, 'none', ('none', 'total', 'atomic', 'total direction', 'atomic direction'), None],
            'fixed_magnetization': [0, float, 0.0, None, None], #check i (1 .. 3)
            'lambda': [0, float, 1.0, None, None],
            'report': [0, int, 1, isPositive, None],
            'lspinorb': [0, bool, None, None, None], #TODO default not specified in docs
            'assume_isolated': [0, str, 'none', ('none', 'makov-payne', 'martyna-tuckerman', 'esm'), None],
            'esm_bc': [0, str, 'pbc', ('pbc', 'bc1', 'bc2', 'bc3'), None],
            'esm_w': [0, float, 0.0, None, None],
            'esm_efield': [0, float , 0.0, None, None],
            'esm_nfit': [0, int, 4, None, None],
            'vdw_corr': [0, str, 'none', ('grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d', 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler', 'XDM', 'xdm'), None],
            'london': [0, bool, False, None, None],
            'london_s6': [0, float, 0.75, None, None],
            'london_rcut': [0, float, 200.0, None, None],
            'xdm': [0, bool, False, None, None],
            'xdm_a1': [0, float, 0.6836, None, None],
            'xdm_a2': [0, float, 1.5045, None, None],
            'space_group': [0, int, 0, self._rangeSpaceGroup, self._checkSpaceGroup],
            'uniqueb': [0, bool, False, None, None],
            'origin_choice': [0, int, 1, (1, 2, 3), None], # TODO Maybe there are more possible origins?
            'rhombohedral': [0, bool, True, None, None]
        }
        super().__init__(name, keys)
