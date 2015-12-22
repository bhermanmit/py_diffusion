"""Materials module.
"""

import numpy as np


from diffusion.material import Material


class Materials(object):

    """Materials class contains list of materials.
    """

    def __init__(self, num_materials):

        """Materials constructor.
        """

        self._num_materials = num_materials

        self._materials = np.empty((num_materials,), dtype=Material)

    def __getitem__(self, index):

        if index >= self._num_materials:
            raise IndexError("Material index out-of-bounds.")

        return self._materials[index]

    def set_material(self, material, index):

        """Sets a material in the material array.
        """

        if index >= self._num_materials:
            raise IndexError("Material index out-of-bounds.")

        if not isinstance(material, Material):
            raise TypeError("Material is not of correct type.")

        self._materials[index] = material
