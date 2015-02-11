"""Card: Cell Parameters

"""
class CellParameters:
    """
    Card: CELL_PARAMETERS
    Defines the vectors for lattice if ibrav = 0. This Card is optional.

    bohr     : latticle vectors in bhor radii
               alat = sqrt(vec1*vec1)

    angstrom : lattice vectors in angstroms
               alat = sqrt(vec1*vec1)

    alat     : lattice vectors in units of the lattice parameter
               either celldm(1) or a.

    DEFAULT: alat
    """
    options = ["alat", "bhor", "angstrom"]

    def __init__(self):
        self.name = "CELL_PARAMETERS"
        self.option = None
        self.lattice_vec = None

    def add_lattice_vec(self, vec1, vec2, vec3, option="alat"):
        """
        Sets the lattice vectors [vec1, vec2, vec3].

        option - specifies the units to use:
         - bohr
         - angstrom
         - alat
        """
        if not len(vec1) == len(vec2) == len(vec3):
            raise Exception("Must be vector len 3")

        if option not in CellParameters.options:
            error_str = "CELL_PARAMETERS {0} not valid option".format(self.option)
            raise Exception(error_str)

        self.option = option
        self.lattice_vec = [vec1, vec2, vec3]

    def validate(self):
        """
        Validate that class K_POINTS is properly setup
        (to the best of my knowledge)
        """
        if self.option not in CellParameters.options and not (self.option == None and self.lattice_vec == None):
            error_str = "CELL_PARAMETER {0} not valid (should never happen)".format(self.option)
            raise Exception(error_str)

    def __str__(self):
        cellparameter_str = "{0} ({1})\n".format(self.name, self.option)
        v1_str = " ".join(map(str, self.lattice_vec[0]))
        cellparameter_str += QE_TAB + v1_str + "\n"
        v2_str = " ".join(map(str, self.lattice_vec[1]))
        cellparameter_str += QE_TAB + v2_str + "\n"
        v3_str = " ".join(map(str, self.lattice_vec[2]))
        cellparameter_str += QE_TAB + v3_str + "\n"

        return cellparameter_str

