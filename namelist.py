"""
Namelist

These are dictionaries with configuration varaibles.

# TODO implement to check for key "blank" only used when "key" set to true
"""
import os
import re
from collections import Callable

class Namelist:
    """
    Defines the partial class implementation of each namelist

    The keys info dictionary:
    key-name :
      narg = number of args (eg. for ibrav narg=0, for celldm(i) narg=1, etc.)
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
        # We ignore key_index even if it has values
        key, index = self.parseKey(key)

        if not self.isDefinedKey(key):
            error_str = "{0} key: {1} not valid key".format(self.name, key)
            raise Exception(error_str)

        key_narg, key_type, key_default, key_range, key_config, key_doc = self.getKeyInfo(key)

        ## Key name
        keydesc_str = "Key: '{0}'\n".format(key)

        ## Key Number of indicies
        keydesc_str += "Number Indicies\n   {0}".format(key_narg)

        ## Key set value (set by user)
        keydesc_str += "Set Value\n   '{0}'".format(self.getSetKeyValue(key))

        ## Key value type allowed [0]
        keydesc_str += "Type\n   '{0}'".format(key_type.__name__)

        ## Key default value [1]
        keydesc_str += "Default\n"
        if isinstance(key_default, Callable):
            keydesc_str += "   Function Name '{0}'".format(key_default.__name__)
        keydesc_str += "   Value '{0}'\n".format(self.getDefaultKeyValue(key))

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
        if key_config:
            keydesc_str += "   Function Name '{0}'\n".format(key_config.__name__)
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

    def getSetKeyValue(self, key, index=()):
        """Returns the key value if it has been set by the user. Otherwise
        returns None.

        Assumes: valid key

        """
        values = self.keypairs.get(key, None)
        if values and index in values:
            return values[index]
        # This returns either None, or the whole dictionary
        return values

    def getCurrentKeyValue(self, key, index=()):
        """Returns: key value if set by the user. Otherwise the current
        default value of key (getDefaultKeyValue).

        Assumes: valid key

        """
        value = self.getSetKeyValues(key, index)
        if value:
            return value
        else:
            return self.getDefaultKeyValue(key)

    def getDefaultKeyValue(self, key):
        """Returns: current default value of key.

        Assumes: valid key

        **Remember for some keys, their default value changes with
        respect to other keys, or environment variables
        (eg. key:'calculation')**

        """
        value_default = self.keys.get(key)[2]
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
            if len(key_info) != 6:
                error_str = "key '{0}' length wrong"
                raise Exception(error_str.format(key))

            key_narg, key_type, key_default, key_range, key_config, key_doc = key_info

            # Key Type. [str, float, int, bool]
            if key_type not in [str, float, int, bool]:
                error_str = "key '{0}' [0] type wrong"
                raise Exception(error_str.format(key))

            # Key Default Value Type [Key Type, None, or function]
            if not isinstance(key_default, (key_type, Callable, type(None))):
                error_str = "key '{0}' [1] default wrong"
                raise Exception(error_str.format(key))

            # Key Range Type [Tuple list ,None(everything valid), Function]
            if not isinstance(key_range, (tuple, type(None), Callable)):
                error_str = "key '{0}' [2] range wrong"
                raise Exception(error_str.format(key))

            # Key Config Type [None or Function]
            if not isinstance(key_config, (type(None), Callable)):
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
        key_config = self.getKeyInfo(key)[4]
        if key_config:
            status, message = key_config(qe)
            if not status:
                raise Exception(message)

    def parseKey(self, unparsed_key):
        """Keys have many forms:
        <key_name>(int, int, ... )

        """
        # Remove all whitespace charaters (proud of this one :) )
        key = ''.join(_ for _ in unparsed_key if _ not in " \n\r\t\v")
        key_regex = re.compile("^(\w+)$|^(\w+)\((\d+(?:,\d+)*)\)$")
        key_match = key_regex.match(key)

        if not key_match:
            error_str = "malformed key string '{0}'"
            raise Exception(error_str.format(key))

        key_match = key_match.groups()

        if key_match[0]:
            return [key_match[0], ()]
        else:
            return [key_match[1], tuple(map(int, re.findall("\d+", key_match[2])))]


    def validateKeypair(self, keypair):
        """Determines if keypair is in domain and range of allowed
        values. (DOES NOT DETERMINE if the keypair makes sense and is
        properly configured. checkKeypair does this functionaliy).

        Returns: if keypair is in domain and range
                 of namelist allowed values. Exception otherwise.

        Ensures: keypair is in domain and range of given namelist

        """
        key, index, value = keypair
        
        if not len(keypair) == 3 or \
           not isinstance(key, str) or \
           not isinstance(index, tuple):
            error_str = "keypair consist of ['str', value]"
            raise Exception(error_str)

        # Check if key is valid
        if not self.isDefinedKey(key):
            error_str = "{0} key: '{1}' invalid name".format(self.name, key)
            raise Exception(error_str)

        key_narg, key_type, key_default, key_range, key_config, key_doc = self.getKeyInfo(key)

        # Check if key came with correct number of args
        key_narg = self.getKeyInfo(key)[0]
        if len(index) != key_narg:
            error_str = "{0} key: '{1}' expects {2} indicies {3} given"
            raise Exception(error_str.format(self.name, key, key_narg, index))

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
        # We have to treat indicies differently
        for key, values in self.keypairs.items():
            for index, value in values.items():
                self.validateKeypair([key, index, value])

        # Check that for each key in namelist it is properly
        # configured with respect to its global environment.  (can
        # fail often due to bad user configuration) But that is what
        # this is for anyways!!
        for key in self.keys:
            self.validateKeyConfig(key, qe)


    def addKeypairs(self, keypairs):
        """Adds List, Tuple, or Dictionary of unformated keypairs to the
        namelist.  Exception otherwise. Keypairs up to Exception are
        still added.

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
        key_unparsed, value = keypair
        key, index = self.parseKey(key_unparsed)

        self.validateKeypair([key, index, value])

        # Check if user is setting key to its default value (harmless)
        if value == self.getDefaultKeyValue(key):
            print("Warning: Setting Key '{0}' to its default value".format(key))

        if self.getSetKeyValue(key):
            # Check if user is overwritting key (harmless)
            values = self.getSetKeyValue(key)
            if index in values:
                print("Warning: Overwritting Key '{0}'".format(key))
            values.update({index: value})
        else:
            # Not values currently in keypair
            # This will start a dictionary as well
            self.keypairs[key] = {index: value}

    def keypairToString(self, keypair):
        """Assumes: keypair is a valid Keypair
        of form:

        key, index, value = keypair
        key(index) = value

        Returns: namelist string representation of keypair

        """
        key, index, value = keypair

        keypair_str = "{0}".format(key)
        if index:
            keypair_str += "({0})".format(",".join(map(str, index)))
        keypair_str += " = "
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
        for key, values in self.keypairs.items():
            for index, value in values.items():
                namelist_str += self.NAMELIST_SPACE
                namelist_str += self.keypairToString([key, index, value])
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
        if self.getCurrentKeyValue('calculation') in ['scf', 'nscf', None]:
            return 1
        else:
            return 50

    def _defaultTprnfor(self):
        return self.getCurrentKeyValue('calculation') in ['vc-md', 'vc-relax']

    def _defaultOutdir(self):
        return os.environ.get('ESPRESSO_TMPDIR', './')

    def _defaultPseudoDir(self):
        return os.environ.get('ESPRESSO_PSEUDODIR',
                              os.environ.get('HOME') + '/espresso/pseudo/')

    def _defaultDiskIO(self):
        # Include None because scf is default value
        if self.getCurrentKeyValue('calculation') in ['scf', None]:
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
            'starting_ns_eigenvalue': [0, float, -1.0, None, None], #TODO valid range TODO check
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

        super().__init__(name, keypairs, keys)
