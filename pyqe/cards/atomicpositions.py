"""
Card: Atomic Positions

"""

class AtomicPositions:
    """
    ATOMIC_POSITIONS Card

    options:
    alat       : atomic positions are in cartesian coordinates, in
                 units of the lattice parameter (either celldm(1)
                 or A). If no option is specified, 'alat' is assumed;
                 not specifying units is DEPRECATED and will no
                 longer be allowed in the future

    bohr       : atomic positions are in cartesian coordinate,
                 in atomic units (i.e. Bohr radii)

    angstrom   : atomic positions are in cartesian coordinates,
                 in Angstrom

    crystal    : atomic positions are in crystal coordinates, i.e.
                 in relative coordinates of the primitive lattice
                 vectors as defined either in card CELL_PARAMETERS
                 or via the ibrav + celldm / a,b,c... variables

    crystal_sg : atomic positions are in crystal coordinates, i.e.
                 in relative coordinates of the primitive lattice.
                 This option differs from the previous one because
                 in this case only the symmetry inequivalent atoms
                 are given. The variable space_group must indicate
                 the space group number used to find the symmetry
                 equivalent atoms. The other variables that control
                 this option are uniqueb, origin_choice, and
                 rhombohedral.
    """
    options = ("alat", "bhor", "crystal", "crystal_sg", "angstrom")

    def __init__(self):
        self.name = "ATOMIC_POSITIONS"
        self.option = None
        self.atom_positions = []

    def add_atom_position(self, symbol, position, option="alat"):
        """
        Add atom position:
        symbol    (1-2) characters
        position  [x, y, z]
        """
        if len(symbol) > 2:
            error_str = "Chemical Symbol 1-2 Characters [{0}]".format(symbol)
            raise Exception(error_str)

        if not len(position) == 3:
            raise Exception("Must be vector len 3")

        if option not in AtomicPositions.options:
            error_str = "CELL_PARAMETERS {0} not valid option".format(self.option)
            raise Exception(error_str)

        self.option = option
        self.atom_positions.append([symbol, position])

    def validate(self):
        """
        Validate the class ATOMIC_POSITIONS is properly setup
        (to the best of my knowledge)
        """
        if self.option not in AtomicPositions.options:
            error_str = "ATOMIC_POSITIONS {0} not valid (should never happen)".format(self.option)
            raise Exception(error_str)

    def __str__(self):
        import pyqe.config as config

        atomicpositions_str = "{0} ({1})\n".format(self.name, self.option)
        for atom_position in self.atom_positions:
            symbol = atom_position[0]
            position_str = " ".join(map(str, atom_position[1]))
            atomicpositions_str += config.card_space + symbol + " " + position_str + "\n"
        return atomicpositions_str
