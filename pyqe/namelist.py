"""
Namelist

These are dictionaries with configuration varaibles.

# TODO implement to check for key "blank" only used when "key" set to true
"""
import re
from collections import Callable

from pyqe.docs.pwdocs import addDocsToKeys

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
            error_str = "{0} key: '{1}' value '{2}' not type {3}"
            raise Exception(error_str.format(self.name, key, value, key_type.__name__))

        # Check the range of value is correct
        # Formats:
        # () -> all ranges accepted
        # func(value) -> return True/False whether in range
        # (i, j, ...) -> tuple of accepted values
        if isinstance(key_range, Callable):
            if not key_range(value):
                error_str = "'{0}' key '{1}' value '{2}' invalid range"
                raise Exception(error_str.format(self.name, key, value))
        elif key_range and value not in key_range:
            error_str = "'{0}' key '{1}' value '{2}' invalid range"
            raise Exception(error_str.format(self.name, key, value))


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
