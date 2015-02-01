"""
Namelist

These are dictionaries with configuration varaibles.

# TODO implement to check for key "blank" only used when "key" set to true
"""
import os
from collections import Callable

class Namelist:
    """
    Defines the partial class implementation of each namelist

    The keys info dictionary:
    key-name :
      _type = [ str, float, int, bool ]
      default value = [ function(self) returns value of _type, value of _type, None ]
      range = [ (), (value1, value2, ...), function(value) ]
      config = [ function(global QE object) returns [True/False, str] , None | doc string ]
      doc string = str describing key
    """

    NAMELIST_SPACE = "    "

    def __init__(self, name, keypairs, keys):
        self.name = name
        self.keypairs = keypairs

        from pwdocs import addDocsToKeys
        addDocsToKeys(name, keys)
        self.keys = keys
        self.validateKeyDescriptions()

    def getKeyDescription(self, key):
        """Returns: a nicely formated string description of the key for
        Espresso calculations. Otherwise raises Exception

        Ensures: key is valid
        """
        if not self.isDefinedKey(key):
            error_str += "{0} key: {1} not valid key".format(self.name, key)
            raise Exception(error_str)

        key_type, key_default, key_range, key_config, key_doc = key_info.getKeyInfo()

        ## Key name
        keydesc_str = "Key: '{0}'\n".format(key)

        ## Key set value (set by user)
        keydesc_str += "Set Value\n   '{0}'".format(self.getSetKeyValue(key))

        ## Key value type allowed [0]
        keydesc_str += "Type\n   '{0}'".format(key_type.__name__)

        ## Key default value [1]
        keydesc_str += "Default\n"
        if isinstance(key_default, Callable):
            keydesc_str += "   Function Name '{0}'".format(key_default.__name__)
        keydesc_str += "   Value '{0}'\n".format(self.getKeyDefaultValue(key))

        ## Key valid range [2]
        keydesc_str += "Valid Range\n"
        if isinstance(key_range, Callable):
            keydesc_str += "   Function Name '{0}'\n".format(key_range.__name__)
        elif key_range:
            keydesc_str += "   Range: [{0}]\n".format(", ".join(map(str, key_range)))
        else:
            keydesc_str += "   Full Range\n"

        ## Key config checker (global) [3]
        keydesc_str += "Config:\n"
        if key_correct:
            keydesc_str += "   Function Name '{0}'\n".format(key_check.__name__)
        else:
            keydesc_str += "   No Config Function (Assumed Valid)\n"

        ## Key 'docstring' [4]
        keydesc_str += "Description:\n" + key_doc

        return keydesc_str

    def getKeyInfo(self, key):
        """Returns: Info regarding key.

        Assumes: valid key

        """
        return self.keys.get(key)

    def getSetKeyValue(self, key):
        """Returns the key value if it has been set by the user. Otherwise
        returns None.

        Assumes: valid key

        """
        return self.keypairs.get(key)

    def getCurrentKeyValue(self, key):
        """Returns: key value if set by the user. Otherwise the current
        default value of key (getDefaultKeyValue).

        Assumes: valid key

        """
        return self.keypairs.get(key, self.getDefaultKeyValue(key))

    def getDefaultKeyValue(self, key):
        """Returns: current default value of key.

        Assumes: valid key

        **Remember for some keys, their default value changes with
        respect to other keys, or environment variables
        (eg. key:'calculation')**

        """
        value_default = self.keys.get(key)[1]
        if isinstance(value_default, Callable):
            return value_default()
        return value_default

    def isDefinedKey(self, key):
        """Returns: if key is defined for given namelist

        """
        if self.keys.get(key):
            return True
        else:
            return False

    def validateKeyDescriptions(self):
        """**For Developer** 
        Checking the of key description arrays. This is
        really only usefull to the developer and ensure that errors
        were not made.  probably should create a test folder that runs
        this rather that have this function.

        """
        for key, key_info in self.keys.items():
            if len(key_info) != 5:
                error_str = "key '{0}' length wrong"
                raise Exception(error_str.format(key))

            key_type, key_default, key_range, key_config, key_doc = key_info

            # Key Type. [str, float, int, bool]
            if key_type not in [str, float, int, bool]:
                error_str = "key '{0}' [0] type wrong"
                raise Exception(error_str.format(key))

            # Key Default Value Type [Key Type, None, or function]
            if not isinstance(key_default, (key_type, Callable, type(None))):
                error_str = "key '{0}' [1] default wrong"
                raise Exception(error_str.format(key))

            # Key Range Type [Tuple list or Function]
            if not isinstance(key_range, (tuple, Callable)):
                error_str = "key '{0}' [2] range wrong"
                raise Exception(error_str.format(key))

            # Key Config Type [None or Function]
            if not isinstance(key_config, (Callable, type(None))):
                error_str = "key '{0}' [3] config wrong"
                raise Exception(error_str.format(key))

            # Key Doc String Type [str]
            if not isinstance(key_doc, str):
                error_str = "key '{0}' [4] doc wrong"
                raise Exception(error_str.format(key))

    def validateKeyConfig(self, key, qe):
        """Determines if key is properly configured. Notice we need the
        global state of the QE configuration. Some keys (eg. 'ibrav')
        require global knowledge. Sometimes the keypair is the default keypair.

        Returns: if keypair is properly with respect to global configuration.

        Ensures: keypair is correct given the global config
        """
        key_config = self.getKeyInfo(key)[3]
        if key_config:
            status, message = key_config(qe)
            if not status:
                raise Exception(message)


    def validateKeypair(self, keypair):
        """Determines if keypair is in domain and range of allowed
        values. (DOES NOT DETERMINE if the keypair makes sense and is
        properly configured. checkKeypair does this functionaliy).

        Returns: if keypair is in domain and range
                 of namelist allowed values. Exception otherwise.

        Ensures: keypair is in domain and range of given namelist

        """
        if not len(keypair) == 2 or not isinstance(keypair[0], str):
            error_str = "keypair consist of ['str', value]"
            raise Exception(error_str)

        key, value = keypair

        # Check if key is valid
        if not self.isDefinedKey(key):
            error_str = "{0} key: '{1}' invalid name".format(self.name, key)
            raise Exception(error_str)

        key_type, key_default, key_range, key_config, key_doc = self.getKeyInfo(key)

        # Check if value is of correct type
        if not isinstance(value, key_type):
            error_str = "{0} key: '{1}' value '{2}' not type {3}".format(self.name, key, value, key_type.__name__)
            raise Exception(error_str)

        # Check the range of value is correct
        if (isinstance(key_range, Callable) and not key_range(value)) and \
           (key_range and value not in key_range):
            error_str = "'{0}' key '{1}' value '{2}' invalid range".format(self.name, key, value)
            raise Exception(error_str)

    def validate(self, qe):
        """Validate keypairs completely.

        Ensures: Namelist is setup properly with respect to global
                 configuration. (You should be able to run after this)
        """
        # This validates the domain and range of each user set keypair
        # (should really not fail)
        for keypair in self.keypairs.items():
            self.validateKeypair(keypair)

        # Check that for each key in namelist it is properly
        # configured with respect to its global environment.  (can
        # fail often due to bad user configuration) But that is what
        # this is for anyways!!
        for key in self.keys:
            self.validateKeyConfig(key, qe)


    def addKeypairs(self, keypairs):
        """Adds List, Tuple, or Dictionary of keypairs to the namelist.
        Exception otherwise. Keypairs up to Exception are still added.

        Ensures: keypairs are valid keypairs are added to dictionary

        """
        if isinstance(keypairs, (list, tuple)):
            for keypair in keypairs:
                self.addKeypair(keypair)
        elif isinstance(keypairs, dict):
            for keypair in keypairs.items():
                self.addKeypair(keypair)
        else:
            raise Exception("Keypairs must be type List, Tuple, or Dictionary")

    def addKeypair(self, keypair):
        """Adds valid List or Tuple keypair to the namelist.  Exception
        otherwise.

        Ensures: keypair is valid
                 keypair is added to dictionary

        """
        self.validateKeypair(keypair)

        key, value = keypair

        # Check if user is setting key to its default value (harmless)
        if value == self.getDefaultKeyValue(key):
            print("Warning: Setting Key '{0}' to its default value".format(key))

        # Check if user is overwritting key (harmless)
        if self.getSetKeyValue(key):
            print("Warning: Overwritting Key '{0}'".format(key))

        self.keypairs[key] = value

    def keypairToString(self, keypair):
        """Assumes: keypair is a valid Keypair

        Returns: namelist string representation of keypair

        """
        key, value = keypair

        keypair_str = "{0} = ".format(key)
        if isinstance(value, str):
            keypair_str += "'{1}'".format(key, value)
        else:
            keypair_str += "{1}".format(key, value)

        return keypair_str

    def __str__(self):
        if len(self.keypairs) == 0:
            return ""

        namelist_str = " &{0}\n".format(self.name)

        # Iterate through key value
        for keypair in self.keypairs.items():
            namelist_str += self.NAMELIST_SPACE
            namelist_str += self.keypairToString(keypair)
            namelist_str += "\n"

        namelist_str += " /\n"

        return namelist_str

# Functions to evaluate range
def isPositive(value):
    return value > 0.0

def isWithinOneOfZero(value):
    return abs(value) <= 1.0

def isBtwZeroOne(value):
    return value >= 0.0 and value <= 1.0

def isGTFour(value):
    return value > 4

class Control(Namelist):
    """Control Namelist.


    """

    def _defaultNStep(self):
        # Include None because scf is default value
        if self.keypairs.get('calculation') in ['scf', 'nscf', None]:
            return 1
        else:
            return 50

    def _defaultTprnfor(self):
        return self.keypairs.get('calculation') in ['vc-md', 'vc-relax']

    def _defaultOutdir(self):
        return os.environ.get('ESPRESSO_TMPDIR', './')

    def _defaultPseudoDir(self):
        return os.environ.get('ESPRESSO_PSEUDODIR',
                              os.environ.get('HOME') + '/espresso/pseudo/')

    def _defaultDiskIO(self):
        # Include None because scf is default value
        if self.keypairs.get('calculation') in ['scf', None]:
            return 'low'
        else:
            return 'medium'

    def _checkDirectory(self, key):
        # Check the set key value by user
        directory = self.getSetKeyValue(key)
        if directory:
            if os.path.isdir(directory):
                return [True, None]
            else:
                error_str = "key {0} -> '{1}' user set directory does not exist"
                return [False, error_str.format(key, directory)]

        # key value was not set by user so check default value
        directory = self.getDefaultKeyValue(key)
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
            'calculation': [str, 'scf', ('scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax'), None],
            'title': [str, '', (), None],
            'verbosity': [str, 'low', ('low', 'high'), None],
            'restart_mode': [ str, 'from_scratch', ('from_scratch', 'restart'), None],
            'wf_collect': [bool, False, (), None],
            'nstep': [int, self._defaultNStep, (), None],
            'iprint': [int, None, (), None],
            'tstress': [bool, False, (), None],
            'tprnfor': [bool, self._defaultTprnfor , (), None],
            'dt': [float, 20.0, isPositive, None],
            'outdir': [str, self._defaultOutdir, (), self._checkOutdir],
            'wfcdir': [str, self._defaultOutdir, (), self._checkWfcdir],
            'prefix': [str, 'pwscf', (), None],
            'lkpoint_dir': [bool, True, (), None],
            'max_seconds': [float, 1.0e7, isPositive, None],
            'etot_conv_thr': [float, 1.0e-4, isPositive, None],
            'forc_conv_thr': [float, 1.0e-4, isPositive, None],
            'disk_io': [str, self._defaultDiskIO, ('none', 'low', 'medium', 'high'), None],
            'pseudo_dir': [str, self._defaultPseudoDir, (), self._checkPseudoDir],
            'tefield': [bool, False, (), None],
            'dipfield': [bool, False, (), None],
            'lelfield': [bool, False, (), None],
            'nberrycyc': [int, 1, isPositive, None],
            'lorbm': [bool, False, (), None],
            'lberry': [bool, False, (), None],
            'gdir': [int, None, (1, 2, 3), None],
            'nppstr': [int, None, isPositive, None]
        }
        super().__init__(name, keypairs, keys)

class System(Namelist):
    """System Namelist

    """

    def _rangeIbrav(self, value):
        return value in ([i for i in range(15)] + [-5])

    def _rangeSpaceGroup(self, value):
        return value in (i for i in range(231))

    def _defaultEcutrho(self):
        return 4.0 * self.getCurrentKeyValue('ecutwfc')

    def _defaultEcutfock(self):
        return self.getCurrentKeyValue('ecutrho')

    def _defaultnqx1(self):
        return self.getCurrentKeyValue('nr1')

    def _defaultnqx2(self):
        return self.getCurrentKeyValue('nr2')

    def _defaultnqx3(self):
        return self.getCurrentKeyValue('nr3')

    def _checkIbrav(self, qe):
        pass

    def _checkCelldm(self, qe):
        pass

    def _checkA(self, qe):
        pass

    def _checkB(self, qe):
        pass

    def _checkC(self, qe):
        pass

    def _checkCosAB(self, qe):
        pass

    def _checkCosAC(self, qe):
        pass

    def _checkCosBC(self, qe):
        pass

    def _checkNat(self, qe):
        pass

    def _checkNtyp(self, qe):
        pass

    def _checkNri(self, qe):
        if not qe.control.getSetKeyValue('nr1') and \
           qe.control.getSetKeyValue('nr2') and \
            qe.control.getSetKeyValue('nr3'):
            return [False, "Must set 'nr1', 'nr2', and 'nr3'"]
        return [True, None]

    def _checkNris(self, qe):
        if not qe.control.getSetKeyValue('nr1s') and \
           qe.control.getSetKeyValue('nr2s') and \
            qe.control.getSetKeyValue('nr3s'):
            return [False, "Must set 'nr1s', 'nr2s', and 'nr3s'"]
        return [True, None]

    def _checkSpaceGroup(self, qe):
        pass

    def __init__(self):
        name = "SYSTEM"
        keypairs = {}
        keys = {
            'ibrav': [int, None, self._rangeIbrav, self._checkIbrav],
            'celldm(i)': [float, None, isPositive, self._checkCelldm],
            'A': [float, None, isPositive, self._checkA],
            'B': [float, None, isPositive, self._checkB],
            'C': [float, None, isPositive, self._checkC],
            'cosAB': [float, None, isWithinOneOfZero, self._checkCosAB],
            'cosAC': [float, None, isWithinOneOfZero, self._checkCosAC],
            'cosBC': [float, None, isWithinOneOfZero, self._checkCosBC],
            'nat': [int, None, isPositive, self._checkNat],
            'ntyp': [int, None, isPositive, self._checkNtyp],
            'nbnd': [int, None, isGTFour, None], #TODO default nbnd (insulator)
            'tot_charge': [float, 0.0, (), None],
            'tot_magnetization': [float, -1.0, isWithinOneOfZero, None],
            'starting_magnetization(i)': [float, None, isWithinOneOfZero, None],
            'ecutwfc': [float, None, isPositive, None],
            'ecutrho': [float, self._defaultEcutrho, isPositive, None],
            'ecutfock': [float, self._defaultEcutfock, isPositive, None],
            'nr1': [int, None, isPositive, self._checkNri],
            'nr2': [int, None, isPositive, self._checkNri],
            'nr3': [int, None, isPositive, self._checkNri],
            'nr1s': [int, None, isPositive, self._checkNris],
            'nr2s': [int, None, isPositive, self._checkNris],
            'nr3s': [int, None, isPositive, self._checkNris],
            'nosym': [bool, False, (), None],
            'nosym_evc': [bool, False, (), None],
            'noinv': [bool, False, (), None],
            'no_t_rev': [bool, False, (), None],
            'force_symmorphic': [bool, False, (), None],
            'use_all_frac': [bool, False, (), None],
            'occupations': [str, None, ('smearing', 'tetrahedra', 'fixed', 'from_input'), None],
            'one_atom_occupations': [bool, False, (), None],
            'starting_spin_angle': [bool, False, (), None],
            'degauss': [float, 0.0, isPositive, None],
            'smearing': [str, 'gaussian', ('gaussian', 'methfessel-paxton', 'm-p', 'mp', 'mazari-vanderbilt', 'cold', 'm-v', 'mv', 'fermi-dirac', 'f-d', 'fd'), None],
            'nspin': [int, 1, (1, 2, 4), None],
            'noncolin': [bool, False, (), None],
            'ecfixed': [float, 0.0, isPositive, None],
            'qcutz': [float, 0.0, isPositive, None],
            'q2sigma': [float, 0.1, isPositive, None],
            'input_dft': [str, None, (), None],
            'exx_fraction': [float, None, isBtwZeroOne, None],
            'screening_parameter': [float, 0.106, (), None],
            'exxdiv_treatment': [str, 'gygi-baldereshi', ('gygi-baldereschi', 'vcut_spherical', 'vcut_ws', 'none'), None],
            'x_gamma_extrapolation': [bool, True, (), None],
            'ecutvcut': [float, 0.0, isPositive, None],
            'nqx1': [int, self._defaultnqx1, isPositive, None],
            'nqx2': [int, self._defaultnqx2, isPositive, None],
            'nqx3': [int, self._defaultnqx3, isPositive, None],
            'lda_plus_u': [bool, False, (), None],
            'lda_plus_u_kind': [int, 0, (0, 1), None],
            'Hubbard_U(i)': [float, 0.0, (), None],  #TODO valid range TODO check
            'Hubbard_J0(i)': [float, 0.0, (), None], #TODO valid range TODO check
            'Hubbard_alpha(i)': [float, 0.0, (), None], #TODO valid range TODO check
            'Hubbard_beta(i)': [float, 0.0, (), None], #TODO valid range TODO check
            'Hubbard_J(i,ityp)': [float, 0.0, (), None], #TODO valid range TODO check
            'starting_ns_eigenvalue(m,ispin,I)': [float, -1.0, (), None], #TODO valid range TODO check
            'U_projection_type': [str, 'atomic', ('atomic', 'ortho-atomic', 'norm-atomic', 'file', 'pseduo'), None],
            'edir': [int, None, (1, 2, 3), None],
            'emaxpos': [float, 0.5, isBtwZeroOne, None],
            'eopreg': [float, 0.1, isBtwZeroOne, None], 
            'eamp': [float, 0.001, (), None], #TODO range
            'angle1(i)': [float, None, (), None], #TODO check (1 .. ntyp)
            'angle2(i)': [float, None, (), None], #TODO check (1 .. ntyp)
            'constrained_magnetization': [str, 'none', ('none', 'total', 'atomic', 'total direction', 'atomic direction'), None],
            'fixed_magnetization(i)': [float, 0.0, (), None], #check i (1 .. 3)
            'lambda': [float, 1.0, (), None],
            'report': [int, 1, isPositive, None],
            'lspinorb': [bool, None, (), None], #TODO default not specified in docs
            'assume_isolated': [str, 'none', ('none', 'makov-payne', 'martyna-tuckerman', 'esm'), None],
            'esm_bc': [str, 'pbc', ('pbc', 'bc1', 'bc2', 'bc3'), None],
            'esm_w': [float, 0.0, (), None],
            'esm_efield': [float , 0.0, (), None],
            'esm_nfit': [int, 4, (), None],
            'vdw_corr': [str, 'none', ('grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d', 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler', 'XDM', 'xdm'), None],
            'london': [bool, False, (), None],
            'london_s6': [float, 0.75, (), None],
            'london_rcut': [float, 200.0, (), None],
            'xdm': [bool, False, (), None],
            'xdm_a1': [float, 0.6836, (), None],
            'xdm_a2': [float, 1.5045, (), None],
            'space_group': [int, 0, self._rangeSpaceGroup, self._checkSpaceGroup],
            'uniqueb': [bool, False, (), None],
            'origin_choice': [int, 1, (1, 2, 3), None], # TODO Maybe there are more possible origins?
            'rhombohedral': [bool, True, (), None]
        }

        super().__init__(name, keypairs, keys)
