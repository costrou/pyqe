"""Card: KPoints

"""
class KPoints:
    """
    KPoint options (see function for description):
     - tpiba
     - tpiba_b
     - tpiba_c
     - crystal
     - crystal_b
     - crystal_c
     - gamma

    DEFAULT: tpiba
    """

    options = ["tpiba", "tpiba_b", "tpiba_c",
               "crystal", "crystal_b", "crystal_c",
               "gamma", "automatic"]

    def __init__(self):
        self.name = "K_POINTS"
        self.option = None
        self.config = None

    def from_monkhorst_pack(self, grid, offset):
        """
        Automatically generated uniform grid of k-points, i.e,
        generates ( nk1, nk2, nk3 ) grid with ( sk1, sk2, sk3 )
        offset.  nk1, nk2, nk3 as in Monkhorst-Pack grids k1, k2, k3
        must be 0 ( no offset ) or 1 ( grid displaced by half a grid
        step in the corresponding direction ) BEWARE: only grids
        having the full symmetry of the crystal work with
        tetrahedra. Some grids with offset may not work.

        FORMAT:
        K_POINTS (automatic)
          nk1 nk2 nk3 sk1 sk2 sk3
        """
        # Validate Input
        if not len(grid) == len(offset) == 3:
            raise Exception("Must be vector len 3")

        self.option = "automatic"
        self.config = [grid, offset]

    def from_list(self, kpoints):
        """read k-points in cartesian coordinates,
        in units of 2 pi/a (default)
        """
        # Validate Input
        if not hasattr(kpoints, "__iter__"):
            raise Exception("kpoints must be list of kpoints")

        self.option = "tpiba"
        self.config = kpoints

    def validate(self):
        """
        Validate that class K_POINTS is properly setup
        (to the best of my knowledge)
        """
        if self.option not in KPoints.options:
            error_str = "K_POINT {0} not valid (should never happen)".format(self.option)
            raise Exception(error_str)

    def __str__(self):
        import pyqe.config as config

        kpoint_str = "{0} ({1})\n".format(self.name, self.option)

        if self.option == "automatic":
            # A little trick to convert list of int to delimited string
            grid_str = " ".join(map(str, self.config[0]))
            offset_str = " ".join(map(str, self.config[1]))
            kpoint_str += config.card_space + grid_str + " " + offset_str + "\n"
        elif self.option in ["tpiba", "tpiba_b", "tpiba_c", "crystal", "crystal_b", "crystal_c"]:
            kpoint_str += config.card_space + "{0}\n".format(len(self.config))
            for kpoint in self.config:
                kpoint_str += config.card_space + "{0} {1} {2} 1.0\n".format(kpoint[0], kpoint[1], kpoint[2])
        elif self.option == "gamma":
            # Do nothing for gamma point calculation
            pass
        else:
            error_str = "K_POINT ({0}) not implemeted yet! (should not happen)".format(self.option)
            raise Exception(error_str)

        return kpoint_str
