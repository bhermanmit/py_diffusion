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

        self._fine_mesh = None

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

    @property
    def num_parts(self):

        """Number of parts in mesh.
        """

        return self._num_parts

    @property
    def left_boundary(self):

        """The left boundary albedo.
        """

        return self._left_boundary

    @left_boundary.setter
    def left_boundary(self, left_boundary):

        """Set the left boundary condition.
        """

        if isinstance(left_boundary, str):

            left_boundary = self._validate_boundary(left_boundary)

        self._left_boundary = left_boundary

    @property
    def right_boundary(self):

        """The right boundary albedo.
        """

        return self._right_boundary

    @right_boundary.setter
    def right_boundary(self, right_boundary):

        """Set the right boundary condition.
        """

        if isinstance(right_boundary, str):

            right_boundary = self._validate_boundary(right_boundary)

        self._right_boundary = right_boundary

    def generate(self):

        """Generates the fine mesh.
        """

        self._fine_mesh = self._part_sizes / self._part_mesh

    def get_material(self, index):

        """Retrieves the material for a mesh index.
        """ 

        coarse_index = self._fine_to_coarse(index)

        return self._material_map[coarse_index]

    def get_width(self, index):

        """Retrieves the width for a mesh index.
        """ 

        coarse_index = self._fine_to_coarse(index)

        return self._fine_mesh[coarse_index]

    def extract_x(self):

        """Calcultes the x mesh positions.
        """

        slab_x = np.zeros(self._part_mesh.sum())

        count = 0
        start = 0.0
        for i in xrange(self._num_parts):
            width = self._fine_mesh[i]
            for j in xrange(self._part_mesh[i]):
                slab_x[count] = start + width/2.0 + float(j)*width
                count += 1
            start += self._part_sizes[i]

        return slab_x

    @staticmethod
    def _validate_boundary(boundary):

        """Validates the boundary condition option.
        """

        if boundary == "reflective":
            boundary = 1.0
        elif boundary == "vacuum":
            boundary = 0.0
        elif boundary == "zero":
            boundary = -1.0
        else:
            raise ValueError("Boundary conditions not valid.")

        return boundary

    def _fine_to_coarse(self, fine_index):

        """Returns the coarse mesh index.
        """

        boundary = 0
        for i in xrange(len(self._part_mesh)):
            boundary += self._part_mesh[i]

            if fine_index < boundary:
                return i

        raise Exception("Could not determine coarse index.")
