"""
Namelist

These are dictionaries with configuration varaibles.

# TODO implement to check for key "blank" only used when "key" set to true
"""
import re
from collections import Callable, defaultdict

from pyqe.docs.pwdocs import getPWDocForKey

class KeyInfo():
    """
    Format for storing information about each key
           key-name :  key = name of the key
           narg = number of args (eg. for ibrav narg=0, for celldm(i) narg=1, etc.)
           type = [ str, float, int, bool ]
           default value = [ function(self) returns value of _type, value of _type, None ]
           range = [ (), (value1, value2, ...), function(value) ]
           config = [ function(global QE object) returns [True/False, str] , None | doc string ]
           doc string = str describing key
    """
    def __init__(self, key, narg, _type, default, _range, config, doc):
        self.key = key
        self.narg = narg
        self.type = _type
        self.default = default
        self.range = _range
        self.config = config
        self.doc = doc

    def validate(self):
        """**For Developer**
        Checking the of key description arrays. This is
        really only usefull to the developer and ensure that errors
        were not made.  probably should create a test folder that runs
        this rather that have this function.

        """
        # Key Type. [str, float, int, bool]
        if self.type not in [str, float, int, bool]:
            error_str = "key '{0}' [0] type wrong"
            raise Exception(error_str.format(keyinfo.key))

        # Key Default Value Type [Key Type, None, or function]
        if not isinstance(self.default, (self.type, Callable, type(None))):
            error_str = "key '{0}' [1] default wrong"
            raise Exception(error_str.format(self.key))

        # Key Range Type [Tuple list ,None(everything valid), Function]
        if not isinstance(self.range, (tuple, type(None), Callable)):
            error_str = "key '{0}' [2] range wrong"
            raise Exception(error_str.format(keyinfo.key))

        # Key Config Type [None or Function]
        if not isinstance(self.config, (type(None), Callable)):
            error_str = "key '{0}' [3] config wrong"
            raise Exception(error_str.format(keyinfo.key))

        # Key Doc String Type [str]
        if not isinstance(self.doc, str):
            error_str = "key '{0}' [4] doc wrong"
            raise Exception(error_str.format(keyinfo.key))


    def to_string(self, namelist, index):
        """Returns: a nicely formated string description of the key for
        Espresso calculations.

        """

        keyinfo_str = (
            "Key               {0}"
            "Number Indicies   {1}"
            "Type              {2}"
            "Set Value         {3}"
            "Default Value"
            "{4}"
            "Valid Range"
            "{5}"
            "Config"
            "{6}"
            "Documentation"
            "{7}"
        )

        ## Key Default
        default_str = "   Value '{0}'\n".format(namelist.get_default_value(key))
        if isinstance(self.default, Callable):
            func_name = self.default.__name__
            default_str += "   Function Name '{0}'\n".format(func_name)

        ## Key Range
        if isinstance(self.range, Callable):
            func_name = self.range.__name__
            range_str = "   Function Name '{0}'\n".format(func_name)
        elif key_range:
            range_str = "   Range: {0}\n".format(keyinfo.range)
        else:
            range_str = "   Full Range\n"

        ## Key Config
        if self.config:
            func_name = self.config.__name__
            config_str = "   Function Name '{0}'\n".format(func_name)
        else:
            config_str = "   No Config Function (Assumed Valid)\n"

        return keyinfo_str.format(self.key,
                                  self.narg,
                                  self.type.__name__,
                                  namelist.get_set_value(self.key, index),
                                  default_str,
                                  range_str,
                                  config_str,
                                  self.doc)

    def __str__(self):
        keyinfo_str = "<key: {0} narg: {1} type: {2} >"
        return keyinfo_str.format(self.key, self.narg, self.type.__name__)


class KeyPair():
    """
    A temporary internal representation of a keypair
    key -> index -> value
    """
    def __init__(self, key, index, value):
        self.key = key
        self.index = index
        self.value = value

    def validate(self, namelist):
        """Determines if keypair is in domain and range of allowed
        values. (DOES NOT DETERMINE if the keypair makes sense and is
        properly configured. checkKeypair does this functionaliy).

        Returns: if keypair is in domain and range
                 of namelist allowed values. Exception otherwise.

        Ensures: keypair is in domain and range of given namelist

        """
        if not isinstance(self.key, str) or \
           not isinstance(self.index, tuple):
            error_str = "keypair consist of [str, tuple, value]"
            raise Exception(error_str)

        keyinfo = namelist.get_key_info(self.key)

        # Check if key is valid
        if not keyinfo:
            error_str = "{0} key: '{1}' invalid name"
            raise Exception(error_str.format(namelist.name, self.key))

        # Check if key came with correct number of args
        if len(self.index) != keyinfo.narg:
            error_str = "{0} key: '{1}' expects {2} indicies {3} given"
            raise Exception(error_str.format(
                namelist.name, self.key, keyinfo.narg, self.index))

        # Check if value is of correct type
        if not isinstance(self.value, keyinfo.type):
            error_str = "{0} key: '{1}' value '{2}' not type {3}"
            raise Exception(error_str.format(
                namelist.name, self.key, self.value, keyinfo.type.__name__))

        # Check the range of value is correct
        # Formats:
        # None -> all ranges accepted
        # func(value) -> return True/False whether in range
        # (i, j, ...) -> tuple of accepted values
        if isinstance(keyinfo.range, Callable):
            if not keyinfo.range(self.value):
                error_str = "'{0}' key '{1}' value '{2}' invalid range"
                raise Exception(error_str.format(
                    namelist.name, self.key, self.value))

        elif keyinfo.range and self.value not in keyinfo.range:
            error_str = "'{0}' key '{1}' value '{2}' invalid range"
            raise Exception(error_str.format(
                namelist.name, self.key, self.value))

    def to_string(self):
        """Assumes: keypair is a valid Keypair
        of form:

        key, index, value = keypair
        key(index) = value

        Returns: namelist string representation of keypair

        """
        keypair_str = "{0}".format(self.key)
        if self.index:
            keypair_str += "({0})".format(','.join(map(str, self.index)))
        keypair_str += " = "
        if isinstance(self.value, str):
            keypair_str += "'{0}'".format(self.value)
        else:
            keypair_str += "{0}".format(self.value)

        return keypair_str


class Namelist:
    """
    Defines the partial class implementation of each namelist

    keys is a dictionary of KeyInfo
    """

    NAMELIST_SPACE = "    "

    def __init__(self, name, keys):
        self.name = name
        self.keypairs = defaultdict(dict)
        self.keys = {}

        for key, [narg, _type, default, _range, config] in keys.items():

            keyinfo = KeyInfo(key,
                              narg,
                              _type,
                              default,
                              _range,
                              config,
                              getPWDocForKey(name, key))
            keyinfo.validate()
            self.keys.update({key: keyinfo})

    def describe_key(self, key):
        key, index = self.parse_key(key)

        if not self.keys.get(key):
            error_str = "{0} key: {1} not valid key"
            raise Exception(error_str.format(self.name, key))

    def get_key_info(self, key):
        """Returns the keyinfo object regarding key.

        """
        return self.keys.get(key)

    def get_set_value(self, key, index=()):
        """Returns the key => index value if it has been set by the
        user. Otherwise returns None.

        Assumes: valid key

        """
        values = self.keypairs.get(key)
        if values and index in values:
            return values[index]
        # This returns either None, or the whole dictionary
        return None

    def get_current_value(self, key, index=()):
        """Returns: key -> index -> value if set by the user. Otherwise the
        current default value of key (get_default_value).

        Assumes: valid key

        """
        value = self.get_set_value(key, index)
        if value:
            return value
        return self.get_default_value(key)

    def get_default_value(self, key):
        """Returns: current default value of key.

        Assumes: valid key

        **Remember for some keys, their default value changes with
        respect to other keys, or environment variables
        (eg. key:'calculation')**

        """
        keyinfo = self.keys.get(key)
        if isinstance(keyinfo.default, Callable):
            return keyinfo.default()
        return keyinfo.default

    def validate_config(self, key, qe):
        """Determines if key is properly configured. Notice we need the
        global state of the QE configuration. Some keys (eg. 'ibrav')
        require global knowledge. Sometimes the keypair is the default keypair.

        Returns: if keypair is properly with respect to global configuration.

        Ensures: keypair is correct given the global config
        """
        keyinfo = self.keys.get(key)
        if keyinfo.config:
            result, message = keyinfo.config(qe)
            if not result:
                raise Exception(message)

    def parse_key(self, unparsed_key):
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
            index = map(int, re.findall("\d+", key_match[2]))
            return [key_match[1], tuple(index)]


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
                keypair = KeyPair(key, index, value)
                keypair.validate(self)

        # Check that for each key in namelist it is properly
        # configured with respect to its global environment.  (can
        # fail often due to bad user configuration) But that is what
        # this is for anyways!!
        for key in self.keys:
            self.validate_config(key, qe)


    def add_keypairs(self, keypairs):
        """Adds List, Tuple, or Dictionary of unformated keypairs to the
        namelist.  Exception otherwise. Keypairs up to Exception are
        still added.

        Ensures: keypairs are valid keypairs are added to dictionary

        """
        if isinstance(keypairs, (list, tuple)):
            for keypair in keypairs:
                self.add_keypair(keypair)
        elif isinstance(keypairs, dict):
            for keypair in keypairs.items():
                self.add_keypair(keypair)
        else:
            raise Exception("Keypairs must be type List, Tuple, or Dictionary")

    def add_keypair(self, keypair):
        """Adds valid List or Tuple keypair to the namelist.  Exception
        otherwise.

        Ensures: keypair is valid
                 keypair is added to dictionary

        """
        key_str, value = keypair
        key, index = self.parse_key(key_str)

        keypair = KeyPair(key, index, value)
        keypair.validate(self)

        # Check if user is setting key to its default value (harmless)
        # if value == self.get_default_value(key):
        #     print("Warning: Setting Key '{0}' to its default value".format(key))

        # Check if user is overwritting key (harmless)
        # if self.get_set_value(key, index):
        #     print("Warning: Overwritting Key '{0}'".format(key))

        self.keypairs[key].update({index: value})


    def to_string(self):
        if len(self.keypairs) == 0:
            return ""

        namelist_str = " &{0}\n".format(self.name)

        # Iterate through key value
        for key, values in self.keypairs.items():
            for index, value in values.items():
                keypair = KeyPair(key, index, value)
                namelist_str += self.NAMELIST_SPACE
                namelist_str += keypair.to_string()
                namelist_str += "\n"
        namelist_str += " /\n"

        return namelist_str

    def __str__(self):
        namelist_str = "<Namelist: {0} KeysPairs: {1} Keys: {2}>"
        return namelist_str.format(
            self.name, len(self.keypairs), len(self.keys))
