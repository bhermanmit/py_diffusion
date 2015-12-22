"""Solver module.
"""


class Solver(object):

    """Solver class.
    """

    def __init__(self, materials, mesh):

        """Solver constructor.
        """

        self._materials = materials
        self._mesh = mesh

        self._loss_matrix = None
        self._production_matrix = None

        self._eigenvalues = None
        self._eigenvectors = None

    def build_loss_matrix(self):

        """Builds left hand side loss matrix.
        """

        pass

    def build_production_matrix(self):

        """Builds right hand side production matrix.
        """

        pass

    def solve(self):

        """Solves the diffusion system.
        """

        pass
