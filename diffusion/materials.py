"""Materials module.
"""

import numpy as np


from diffusion.material import Material


class Materials(object):

    """Materials class contains list of materials.
    """

    def __init__(self, num_materials, num_energy_groups):

        """Materials constructor.
        """

        self._num_materials = num_materials
        self._num_energy_groups = num_energy_groups

        self._materials = np.empty((num_materials,), dtype=Material)

    def __getitem__(self, index):

        if index >= self._num_materials:
            raise IndexError("Material index out-of-bounds.")

        return self._materials[index]

    @property
    def num_materials(self):

        """The number of materials.
        """

        return self._num_materials

    @property
    def num_energy_groups(self):

        """The number of energy groups.
        """

        return self._num_energy_groups

    def set_material(self, material, index):

        """Sets a material in the material array.
        """

        if index >= self._num_materials:
            raise IndexError("Material index out-of-bounds.")

        if not isinstance(material, Material):
            raise TypeError("Material is not of correct type.")

        if self._num_energy_groups != material.num_energy_groups:
            raise ValueError("Material energy groups do not match.")

        self._materials[index] = material

    def generate(self):

        """Generates all the materials.
        """

        for material in self._materials:
            material.generate()
