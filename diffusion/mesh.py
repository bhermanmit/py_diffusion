"""Mesh module.
"""

import numpy as np


class Mesh(object):

    """Mesh class.
    """

    def __init__(self, num_parts):

        """Mesh constructor.
        """

        self._num_parts = num_parts
        self._material_map = None
        self._part_sizes = None
        self._part_mesh = None

        self._left_boundary = ''
        self._right_boundary = ''

    @property
    def material_map(self):

        """Map of material ids for parts.
        """

        return self._material_map

    @material_map.setter
    def material_map(self, material_map):

        """Sets the material map.
        """

        material_map = np.array(material_map)

        if material_map.shape[0] != self._num_parts:
            raise ValueError("Size of material map invalid!")

        self._material_map = material_map

    @property
    def part_sizes(self):

        """Size of each part in cm.
        """

        return self._part_sizes

    @part_sizes.setter
    def part_sizes(self, part_sizes):

        """Sets the part sizes.
        """

        part_sizes = np.array(part_sizes)

        if part_sizes.shape[0] != self._num_parts:
            raise ValueError("Size of part sizes invalid!")

        self._part_sizes = part_sizes

    @property
    def part_mesh(self):

        """The number of mesh bins per part.
        """

        return self._part_mesh

    @part_mesh.setter
    def part_mesh(self, part_mesh):

        """Sets the part mesh.
        """

        part_mesh = np.array(part_mesh)

        if part_mesh.shape[0] != self._num_parts:
            raise ValueError("Size of part mesh invalid!")

        self._part_mesh = part_mesh

    def set_left_boundary(self, left_boundary):

        """Set the left boundary condition.
        """

        self._validate_boundary(left_boundary)

        self._left_boundary = left_boundary

    def set_right_boundary(self, right_boundary):

        """Set the right boundary condition.
        """

        self._validate_boundary(right_boundary)

        self._right_boundary = right_boundary

    @staticmethod
    def _validate_boundary(boundary):

        """Validates the boundary condition option.
        """

        if boundary != "reflective" and boundary != "vacuum" and \
            boundary != "zero":

            raise ValueError("Boundary conditions not valid.")
