"""
Card Atomic Species

"""
import numpy as np

QE_TAB = "   "

class AtomicSpecies:
    """
    Card: ATOMIC_SPECIES

    Atomic Species defines the atoms to be used in the dft calculation.
    For each line the symbol, mass, and pseudo potential to be used for the
    atom is supplied.
    """
    def __init__(self):
        self.name = "ATOMIC_SPECIES"
        self.atoms = []

    def num_atoms(self):
        """
        returns the numbers of atom types defined
        """
        return len(self.atoms)

    def add_atom_type(self, symbol, mass, pseudopot):
        """
        Add atom:
        symbol    (1-2) characters
        mass      double
        pseudopot path to pseudo potential file
        """
        atom_type = [symbol, mass, pseudopot]

        self.validate_atom_type(atom_type)
        self.atoms.append(atom_type)

    def validate_atom_type(self, atom_type):
        symbol, mass, pseudopot = atom_type

        if len(symbol) > 2:
            error_str = "Chemical Symbol 1-2 Characters [{0}]"
            raise Exception(error_str.format(symbol))

        if not isinstance(mass, float):
            error_str = "Mass {0} not of correct type double"
            raise Exception(error_str.format(mass))

        if not isinstance(pseudopot, str):
            error_str = "PsuedoPot filename {0} must be string"
            raise Exception(error_str.format(pseudopot))

    def validate(self):
        """
        Validate atoms
        """
        for atom in self.atoms:
            self.validate_atom_type(atom)

    def __str__(self):
        atoms_str = "{0}\n".format(self.name)
        for atom in self.atoms:
            atoms_str += QE_TAB + "{0} {1} {2}\n".format(atom[0], atom[1], atom[2])
        return atoms_str
