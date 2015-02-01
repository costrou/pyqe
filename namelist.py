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

    def __init__(self):
        self.keys = {}
        self.keypairs = {}
        self.name = ""

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
        self.name = "CONTROL"
        self.keypairs = {}
        self.keys = {
            'calculation': (str, 'scf', ('scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax'), None, """
A string describing the task to be performed
vc = variable cell
"""),
            'title': (str, '', (), None, """
reprinted on output
"""),
            'verbosity': (str, 'low', ('low', 'high'), None, """
Currently two verbosity levels are implemented:
  'high' and 'low'. 'debug' and 'medium' have the same
  effect as 'high'; 'default' and 'minimal', as 'low'
"""),
            'restart_mode': ( str, 'from_scratch', ('from_scratch', 'restart'), None, """
'from_scratch'  : from scratch. This is the normal way
                  to perform a PWscf calculation
'restart'       : from previous interrupted run. Use this
                  switch only if you want to continue an
                  interrupted calculation, not to start a
                  new one, or to perform non-scf calculations.
                  Works only if the calculation was cleanly
                  stopped using variable "max_seconds", or
                  by user request with an "exit file" (i.e.:
                  create a file "prefix".EXIT, in directory
                  "outdir"; see variables "prefix", "outdir")
"""),
            'wf_collect': (bool, False, (), None, """
This flag controls the way wavefunctions are stored to disk :

True    collect wavefunctions from all processors, store them
        into the output data directory "outdir"/"prefix".save,
        one wavefunction per k-point in subdirs K000001/,
        K000001/, etc.. Use this if you want wavefunctions
        to be readable on a different number of processors.

False   do not collect wavefunctions, leave them in temporary
        local files (one per processor). The resulting format
        will be readable only by jobs running on the same
        number of processors and pools. Requires less I/O
        than the previous case.

Note that this flag has no effect on reading, only on writing.
"""),
            'nstep': (int, self._defaultNStep, (), None, """
number of ionic + electronic steps
"""),
            'iprint': (int, None, (), None, """
band energies are written every "iprint" iterations
"""),
            'tstress': (bool, False, (), None, """
calculate stress. It is set to True automatically if
calculation='vc-md' or 'vc-relax'
"""),
            'tprnfor': (bool, self._defaultTprnfor , (), None, """
calculate forces. It is set to .TRUE. automatically if
calculation='relax','md','vc-md'
"""),
            'dt': (float, 20.0, isPositive, None, """
time step for molecular dynamics, in Rydberg atomic units
(1 a.u.=4.8378 * 10^-17 s : beware, the CP code uses
 Hartree atomic units, half that much!!!)
"""),
            'outdir': (str, self._defaultOutdir, (), self._checkOutdir, """
input, temporary, output files are found in this directory,
see also "wfcdir"
"""),
            'wfcdir': (str, self._defaultOutdir, (), self._checkWfcdir, """
this directory specifies where to store files generated by
each processor (*.wfc{N}, *.igk{N}, etc.). Useful for
machines without a parallel file system: set "wfcdir" to
a local file system, while "outdir" should be a parallel
or networkfile system, visible to all processors. Beware:
in order to restart from interrupted runs, or to perform
further calculations using the produced data files, you
may need to copy files to "outdir". Works only for pw.x.
"""),
            'prefix': (str, 'pwscf', (), None, """
prepended to input/output filenames:
prefix.wfc, prefix.rho, etc.
"""),
            'lkpoint_dir': (bool, True, (), None, """
If .false. a subdirectory for each k_point is not opened
in the "prefix".save directory; Kohn-Sham eigenvalues are
stored instead in a single file for all k-points. Currently
doesn't work together with "wf_collect"
"""),
            'max_seconds': (float, 1.0e7, isPositive, None, """
jobs stops after "max_seconds" CPU time. Use this option
in conjunction with option "restart_mode" if you need to
split a job too long to complete into shorter jobs that
fit into your batch queues. Default is 150 days.
"""),
            'etot_conv_thr': (float, 1.0e-4, isPositive, None, """
convergence threshold on total energy (a.u) for ionic
minimization: the convergence criterion is satisfied
when the total energy changes less than "etot_conv_thr"
between two consecutive scf steps. Note that "etot_conv_thr"
is extensive, like the total energy.
See also "forc_conv_thr" - both criteria must be satisfied
"""),
            'forc_conv_thr': (float, 1.0e-4, isPositive, None, """
convergence threshold on forces (a.u) for ionic minimization:
the convergence criterion is satisfied when all components of
all forces are smaller than "forc_conv_thr".
See also "etot_conv_thr" - both criteria must be satisfied
"""),
            'disk_io': (str, self._defaultDiskIO, ('none', 'low', 'medium', 'high'), None, """
Specifies the amount of disk I/O activity
'high':   save all data to disk at each SCF step

'medium': save wavefunctions at each SCF step unless
          there is a single k-point per process (in which
          case the behavior is the same as 'low')

'low' :   store wfc in memory, save only at the end

'none':   do not save anything, not even at the end
          ('scf', 'nscf', 'bands' calculations; some data
           may be written anyway for other calculations)

Default is 'low' for the scf case, 'medium' otherwise.
Note that the needed RAM increases as disk I/O decreases!
It is no longer needed to specify 'high' in order to be able
to restart from an interrupted calculation (see "restart_mode")
but you cannot restart in disk_io='none'
"""),
            'pseudo_dir': (str, self._defaultPseudoDir, (), self._checkPseudoDir, """
directory containing pseudopotential files
"""),
            'tefield': (bool, False, (), None, """
If True a saw-like potential simulating an electric field
is added to the bare ionic potential. See variables "edir",
"eamp", "emaxpos", "eopreg" for the form and size of
the added potential.
"""),
            'dipfield': (bool, False, (), None, """
If True and tefield=True a dipole correction is also
added to the bare ionic potential - implements the recipe
of L. Bengtsson, PRB 59, 12301 (1999). See variables "edir",
"emaxpos", "eopreg" for the form of the correction. Must
be used ONLY in a slab geometry, for surface calculations,
with the discontinuity FALLING IN THE EMPTY SPACE.
"""),
            'lelfield': (bool, False, (), None, """
If .TRUE. a homogeneous finite electric field described
through the modern theory of the polarization is applied.
This is different from "tefield=.true." !
"""),
            'nberrycyc': (int, 1, isPositive, None, """
In the case of a finite electric field  ( lelfield == .TRUE. )
it defines the number of iterations for converging the
wavefunctions in the electric field Hamiltonian, for each
external iteration on the charge density
"""),

            'lorbm': (bool, False, (), None, """
If .TRUE. perform orbital magnetization calculation.
If finite electric field is applied (lelfield=.true.)
only Kubo terms are computed
[for details see New J. Phys. 12, 053032 (2010)].
The type of calculation is 'nscf' and should be performed
on an automatically generated uniform grid of k points.
Works ONLY with norm-conserving pseudopotentials.
"""),
            'lberry': (bool, False, (), None, """
If .TRUE. perform a Berry phase calculation
See the header of PW/src/bp_c_phase.f90 for documentation
"""),
            'gdir': (int, None, (1, 2, 3), None, """
For Berry phase calculation: direction of the k-point
strings in reciprocal space. Allowed values: 1, 2, 3
1=first, 2=second, 3=third reciprocal lattice vector
For calculations with finite electric fields
(lelfield==.true.) "gdir" is the direction of the field
"""),
            'nppstr': (int, None, isPositive, None, """
For Berry phase calculation: number of k-points to be
calculated along each symmetry-reduced string
The same for calculation with finite electric fields
(lelfield=.true.)
""")}

