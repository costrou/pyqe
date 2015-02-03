"""Control Namelist

"""
from pyqe.namelist import Namelist

import os

# Functions to evaluate range
def isPositive(value):
    return value > 0.0


class Control(Namelist):
    """Control Namelist.


    """

    def _defaultNStep(self):
        # Include None because scf is default value
        if self.get_current_value('calculation') in ['scf', 'nscf', None]:
            return 1
        else:
            return 50

    def _defaultTprnfor(self):
        return self.get_current_value('calculation') in ['vc-md', 'vc-relax']

    def _defaultOutdir(self):
        return os.environ.get('ESPRESSO_TMPDIR', './')

    def _defaultPseudoDir(self):
        return os.environ.get('ESPRESSO_PSEUDODIR',
                              os.environ.get('HOME') + '/espresso/pseudo/')

    def _defaultDiskIO(self):
        # Include None because scf is default value
        if self.get_current_value('calculation') in ['scf', None]:
            return 'low'
        else:
            return 'medium'

    def _checkDirectory(self, key):
        # Check the set key value by user
        directory = self.get_set_value(key)
        if directory:
            if os.path.isdir(directory):
                return [True, None]
            else:
                error_str = "key {0} -> '{1}' user set directory does not exist"
                return [False, error_str.format(key, directory)]

        # key value was not set by user so check default value
        directory = self.get_default_value(key)
        if os.path.isdir(directory):
            return [True, None]

        error_str = "key {0} -> '{1}' default directory does not exist"
        return [False, error_str.format(key, directory)]


    def _checkPseudoDir(self, qe):
        return self._checkDirectory("pseudo_dir")

    def _checkOutdir(self, qe):
        return self._checkDirectory("outdir")

    def _checkWfcdir(self, qe):
        return self._checkDirectory("wfcdir")

    def __init__(self):
        name = "CONTROL"
        keypairs = {}
        keys = {
            'calculation': [0, str, 'scf', ('scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax'), None],
            'title': [0, str, '', None, None],
            'verbosity': [0, str, 'low', ('low', 'high'), None],
            'restart_mode': [0, str, 'from_scratch', ('from_scratch', 'restart'), None],
            'wf_collect': [0, bool, False, None, None],
            'nstep': [0, int, self._defaultNStep, None, None],
            'iprint': [0, int, None, None, None],
            'tstress': [0, bool, False, None, None],
            'tprnfor': [0, bool, self._defaultTprnfor , None, None],
            'dt': [0, float, 20.0, isPositive, None],
            'outdir': [0, str, self._defaultOutdir, None, self._checkOutdir],
            'wfcdir': [0, str, self._defaultOutdir, None, self._checkWfcdir],
            'prefix': [0, str, 'pwscf', None, None],
            'lkpoint_dir': [0, bool, True, None, None],
            'max_seconds': [0, float, 1.0e7, isPositive, None],
            'etot_conv_thr': [0, float, 1.0e-4, isPositive, None],
            'forc_conv_thr': [0, float, 1.0e-4, isPositive, None],
            'disk_io': [0, str, self._defaultDiskIO, ('none', 'low', 'medium', 'high'), None],
            'pseudo_dir': [0, str, self._defaultPseudoDir, None, self._checkPseudoDir],
            'tefield': [0, bool, False, None, None],
            'dipfield': [0, bool, False, None, None],
            'lelfield': [0, bool, False, None, None],
            'nberrycyc': [0, int, 1, isPositive, None],
            'lorbm': [0, bool, False, None, None],
            'lberry': [0, bool, False, None, None],
            'gdir': [0, int, None, (1, 2, 3), None],
            'nppstr': [0, int, None, isPositive, None]
        }
        super().__init__(name, keypairs, keys)
