"""
An assistant to QE (Quantum Espresso) via python
currently my aim is to support the fortran syntax

See
http://www.quantum-espresso.org/wp-content/uploads/Doc/pw_user_guide/node8.html
for details on the input format to PW

'!#' - fortran comment characters


"""
from pyqe.cards import AtomicSpecies, AtomicPositions, KPoints, CellParameters
from pyqe.namelists import Control, System, Electrons, Ions, Cell
from pyqe.io import read_out_file
import sys

class QE:
    """
    Quantum Espresso Main Class

    From this class you can:
     - initialize the input
     - create inputfile
     - run pw.x
    """

    def __init__(self, qe_keypairs):
        self.control = Control()
        self.system = System()
        self.electrons = Electrons()
        self.ions = Ions()
        self.cell = Cell()

        self.atomic_species = AtomicSpecies()
        self.atomic_positions = AtomicPositions()
        self.k_points = KPoints()
        self.cell_parameters = CellParameters()
        # Not Implemented (Card not a class)
        #self.occupations = Card("OCCUPATIONS")
        #self.constrains = Card("CONTRAINTS")
        #self.atomicforces = Card("ATOMIC_FORCES")

        self.namelist_asoc = {
            "control": self.control,
            "system": self.system,
            "electrons": self.electrons,
            "ions": self.ions,
            "cell": self.cell
        }

        self.add_keypairs_to_namespace(qe_keypairs)

    def add_keypairs_to_namespace(self, qe_keypairs):
        """
        Adds the respective keys to each namelist
        """
        for name, keypairs in qe_keypairs.items():
            namelist = self.namelist_asoc.get(name.lower())
            if namelist:
                namelist.add_keypairs(keypairs)
            else:
                error_str = "{0} is not valid namelist"
                raise Exception(error_str.format(name))

    def to_string(self, header=True):
        qe_str = ""

        if (header == True):
            qe_str += "! File Autogenerated from Python QE\n"

        ## NameLists
        qe_str += self.control.to_string()
        qe_str += self.system.to_string()
        qe_str += self.electrons.to_string()
        qe_str += self.ions.to_string()
        qe_str += self.cell.to_string()

        ## Cards
        qe_str += str(self.atomic_species)

        # Only needed if calculations is not
        # 'band' or 'nscf'
        if self.control.get_current_value("calculation") not in ["nscf", "bands"]:
            qe_str += str(self.atomic_positions)

        qe_str += str(self.k_points)

        # Only needed if unitcell is not defined
        # By ibrav
        if self.system.get_current_value("ibrav") == 0:
            qe_str += str(self.cell_parameters)

        # Not Implemented
        # qe_str += str(self.occupations)
        # qe_str += str(self.contraints)
        # qe_str += str(self.atomic_forces)
        return qe_str

    def to_file(self, filename, input_format="fortran"):
        """
        Writes QE configuration to <filename>
        in format specified. Currently only supports
        the Fortran style.
        """
        with open(filename, "w") as qefile:
            qefile.write(self.to_string())

    def run(self, infile="", outfile="", errfile=""):
        """
        Runs QE pw.x.

        If stdin, stdout, stderr filenames are not defined
        no file is created for the given input or output.

        If 'in_filename' is defined the program will run from the
        file rather than stdin via '-i'.

        Notice:
        QE will still create the save files in the directory
        specified by 'outfile' in control namelist
        """
        from subprocess import Popen, PIPE

        prefix = []
        postfix = []

        if infile != "":
            self.to_file(infile)

            pw_command = prefix + ["pw.x", '-i', infile] + postfix
            proc = Popen(pw_command, stdout=PIPE, stderr=PIPE)
            pw_output = proc.communicate()
        else:
            pw_input = self.to_string()

            pw_command = prefix + ["pw.x"] + postfix
            proc = Popen(pw_command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            pw_output = proc.communicate(pw_input.encode())

        proc.wait()
        if proc.returncode != 0:
            with open("CRASH", "r") as f:
                print("Quantum Espresso CRASH FILE:\n{0}".format(f.read()),
                       file=sys.stderr)
            raise Exception("pw.x CRASHED")
                
        pw_out = pw_output[0].decode()
        pw_err = pw_output[1].decode()

        if outfile != "":
            with open(outfile, "w") as f:
                f.write(pw_out)

        if errfile != "":
            with open(errfile, "w") as f:
                f.write(pw_err)

        return read_out_file(pw_out)

    def validate(self):
        """
        Each Namelist and Cards will validate its contents.
        Sometimes they will need access to global information.
        (not sure how to handle this yet)
        """
        self.control.validate(self)
        self.system.validate(self)
        self.electrons.validate(self)
        self.ions.validate(self)
        self.cell.validate(self)

        self.atomic_species.validate()
        self.atomic_positions.validate()
        self.k_points.validate()
        self.cell_parameters.validate()
        # Not Implemented
        # self.occupations.validate()
        # self.constrains.validate()
        # self.atomicforces.validate()


