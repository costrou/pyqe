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

    def read_output(self, output_str):
        """
        Routine I am making currently just for my homework
        but this routine will be expanded. (hack currently)

        Retruns a dictionary of interesting values
        Uses regular expressions extensiely
        """
        import re

        # Extract the information from the header
        header_end = output_str.find("Self-consistent Calculation")
        header_str = output_str[:header_end]

        scf_block_begin = output_str.find("Self-consistent Calculation")
        scf_block_end = output_str.find("convergence has been achieved in")
        scf_block_str = output_str[scf_block_begin:scf_block_end]

        bfgs_steps = []

        bfgs_begin = output_str.find("BFGS Geometry Optimization")
        bfgs_end = output_str.find("End of BFGS Geometry Optimization")
        bfgs_str = output_str[bfgs_begin:bfgs_end]

        if bfgs_str:
            bfgs_block_begin = bfgs_str.find("number of scf cycles")
            bfgs_block_end = bfgs_str.find("Writing output data file")

            if bfgs_block_begin == -1:
                final_begin = output_str.find("Final energy")
                final_end = output_str.find("End final coordinates")
                final_str = output_str[final_begin:final_end]
                bfgs_steps.append((scf_block_str, final_str))
            else:
                bfgs_block_str = bfgs_str[bfgs_block_begin:bfgs_block_end]
                bfgs_steps.append((scf_block_str, bfgs_block_str))

                current_pos = bfgs_block_end + 10 #HACK fix!!!!

                while True:
                    scf_block_begin = bfgs_str[current_pos:].find("Self-consistent Calculation")
                    scf_block_end = bfgs_str[current_pos:].find("convergence has been achieved in")
                    scf_block_str = bfgs_str[current_pos:][scf_block_begin:scf_block_end]

                    bfgs_block_begin = bfgs_str[current_pos:].find("number of scf cycles")
                    bfgs_block_end = bfgs_str[current_pos:].find("Writing output data file")

                    if bfgs_block_begin == -1:
                        final_begin = output_str.find("Final enthalpy")
                        final_end = output_str.find("End final coordinates")
                        final_str = output_str[final_begin:final_end]
                        bfgs_steps.append((scf_block_str, final_str))
                        break

                    bfgs_block_str = bfgs_str[current_pos:][bfgs_block_begin:bfgs_block_end]
                    bfgs_steps.append((scf_block_str, bfgs_block_str))

                    current_pos += bfgs_block_end + 10 #HACK fix!!!!
        else:
            bfgs_steps.append((scf_block_str, ""))

        footer_begin = output_str.find("init_run")
        footer_str = output_str[footer_begin:]

        results = {}

        ## Begin parsing the blocked strings
        # Parse the header string for needed values
        header_values = {
            "bravais-lattice index": (r"bravais-lattice index\s+=\s+(\d+)", float),
            "lattice parameter": (r"lattice parameter \(alat\)\s+=\s+(\d+\.\d+) a\.u\.", float),
            "volume": (r"unit-cell volume\s+=\s+(\d+\.\d+) \(a\.u\.\)\^3", float),
            "number atoms/cell": (r"number of atoms/cell\s+=\s+(\d+)", int),
            "number atom types": (r"number of atomic types\s+=\s+(\d+)", int),
            "number electrons": (r"number of electrons\s+=\s+(\d+\.\d+)", float),
            "number of Kohn Sham states": (r"number of Kohn Sham states\s*=\s+(\d+)", int),
            "kinetic-energy cutoff": (r"kinetic-energy cutoff\s+=\s+(\d+\.\d+)\s+Ry", float),
            "charge density cutoff": (r"charge density cutoff\s+=\s+(\d+\.\d+)\s+Ry", float),
            "convergence threshold": (r"convergence threshold\s+=\s+(\d+\.\d+(?:[eE][-+]?[0-9]+)?)", float),
            "mixing beta": (r"mixing beta\s+=\s+(\d+\.\d+)", float),
            "nstep": (r"nstep\s+=\s=(\d+)", int),
            "celldm": (r"celldm\(1\)=\s+(\d+\.\d+)\s+celldm\(2\)=\s+(\d+\.\d+)\s+celldm\(3\)=\s+(\d+\.\d+)\s+celldm\(4\)=\s+(\d+\.\d+)\s+celldm\(5\)=\s+(\d+\.\d+)\s+celldm\(6\)=\s+(\d+\.\d+)", float),
            "crystal axes": (r"a\(1\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+\)\s+a\(2\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+\)\s+a\(3\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+\)\s+",float),
            "reciprocal axes": (r"b\(1\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+\)\s+b\(2\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+\)\s+b\(3\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+\)\s+",float),
            "FFT dimensions": (r"FFT dimensions:\s+\(\s+(\d+),\s+(\d+),\s+(\d+)\)", int),
        }

        header = {}

        # Find matching keypairs
        for key, (regex, _type) in header_values.items():
            match = re.search(regex, header_str)
            if match:
                if len(match.groups()) == 1:
                    header.update({key: _type(match.groups()[0])})
                else:
                    header.update({key: [_type(_) for _ in match.groups()]})

        # Find all kpoints
        kpoint_regex = re.compile(r"k\(\s+(\d+)\)\s+=\s+\(\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\), wk = \s+(\d+\.\d+)")
        matches = kpoint_regex.findall(header_str)
        kpoints = []
        for match in matches:
            kpoints.append([match[0], [match[1], match[2], match[3]], match[4]])
        header.update({"kpoints": kpoints})

        results.update({"header": header})

        # Parse the results of calculation
        total_energy_regex = re.compile("!\s+total energy\s+=\s+([+-]\d+\.\d+) Ry")
        volume_regex = re.compile("new unit-cell volume\s=\s+(\d+\.\d+) a\.u\.\^3")
        lattice_regex = re.compile("CELL_PARAMETERS.*\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)")
        ion_position_regex = re.compile("([A-Z][a-z]?)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)")

        calculation = {}
        iterations = []
        for scf_block, bgfs_block in bfgs_steps:
            iteration = {}

            match = total_energy_regex.search(scf_block)
            iteration.update({"total energy": float(match.group(1))})

            match = volume_regex.search(bgfs_block)
            if match:
                iteration.update({"volume": float(match.group(1))})

            match = lattice_regex.search(bgfs_block)
            lattice = None
            if match:
                lattice = [float(_) for _ in match.groups()]
                iteration.update({"lattice": [lattice[0:3],
                                              lattice[3:6],
                                              lattice[6:9]]})

            match = ion_position_regex.findall(bgfs_block)
            if match:
                iteration.update({"ion positions": [[_[0], float(_[1]), float(_[2]), float(_[3])] for _ in match]})

            iterations.append(iteration)

        calculation.update({"iterations": iterations})
        calculation.update(iterations[-1])

        results.update({"calculation": calculation})

        return results

    def read_charge_density_file(self, inputfile):
        """
        Reads charge-density.dat file generated by a quantum espresso
        scf/relax/vc-relax run
        """
        import re
        import struct
        import numpy as np

        f = open(inputfile, "rb")
        charge_xml = f.read()

        info_str = b'<INFO nr1="(\d+)" nr2="(\d+)" nr3="(\d+)"/>'
        nrz_start_str = '<z\.{0} type="(\w+)" size="(\d+)" kind="(\d+)">\n'
        nrz_end_str = '\n    </z\.{0}>\n'

        match = re.search(info_str, charge_xml).groups()
        nr1, nr2, nr3 = [int(_.decode()) for _ in match]

        data = []
        for nrz in range(nr3):
            match_start = re.search(nrz_start_str.format(nrz+1).encode(),
                                    charge_xml)
            match_end = re.search(nrz_end_str.format(nrz+1).encode(),
                                  charge_xml)

            # Charge Density Format
            # For each nrz
            # 12 bytes header ... is it from a struct?
            # nr1*nr2 doubles
            # 24 bytes footer ... must be from the struct
            start_offset = 12
            end_offset = 24

            start_index = match_start.end() + start_offset
            end_index = match_end.start() - end_offset

            data += struct.unpack("d"*nr1*nr2,
                                  charge_xml[start_index:end_index])

        charge_data = np.array(data, order='F', ndmin=3)
        charge_data.shape = [nr1, nr2, nr3]
        return charge_data

    def read_eigenvalue_file(self, inputfile):
        """
        Reads the output files for kpoints generated by pw.x
        """
        import xml.etree.ElementTree as ET
        import re
        double_regex = re.compile('[-+]?\d+(?:\.\d*)|\d+\.\d+(?:[eE][-+]?\d+)')

        tree = ET.parse(inputfile)
        root = tree.getroot()

        eigenvalues = [float(_) for _ in
                       double_regex.findall(root.find('EIGENVALUES').text)]
        occupations = [float(_) for _ in
                       double_regex.findall(root.find('OCCUPATIONS').text)]

        return {"eigenvalues": eigenvalues,
                "occupations": occupations}

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

        return self.read_output(pw_out)

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


