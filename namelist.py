"""
Namelist

These are dictionaries with configuration varaibles.
"""
import os
from collections import Callable



class Namelist:
    """
    Defines the partial class implementation of each namelist
    """

    NAMELIST_SPACE = "    "

    def __init__(self, name, keypairs, keys):
        self.name = name
        self.keypairs = keypairs

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
        if (isinstance(key_range, Callable) and not key_range(value)) or \
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

def isPositive(self, value):
    return value > 0.0

class Control(Namelist):
    """
    General variables for controlling the run

    The keys tuple:
    key-name | type | default value | available values (may be function) | doc string for key
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

        from pwdocs import addDocsToKeys
        addDocsToKeys(name, keys)

        super().__init__(name, keypairs, keys)
