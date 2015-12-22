"""Material definition.
"""

import numpy as np


class Material(object):

    """Material class.
    """

    def __init__(self, num_energy_groups):

        """Material constructor.
        """

        self._num_energy_groups = num_energy_groups

        self._total = None
        self._scattering = None
        self._nu_fission = None
        self._diffusion = None
        self._chi = None

    @property
    def total(self):

        """The total macroscopic cross section.
        """

        return self._total

    @total.setter
    def total(self, total):

        """Sets the total macroscopic cross section.
        """

        total = np.array(total)

        if total.shape[0] != self._num_energy_groups:
            raise ValueError("Size of total cross section invalid!")

        self._total = np.array(total)

    @property
    def scattering(self):

        """The scattering macroscopic cross section.
        """

        return self._scattering

    @scattering.setter
    def scattering(self, scattering):

        """Sets the scattering macroscopic cross section.
        """

        scattering = np.array(scattering)

        if len(scattering.shape) != 2:
            raise ValueError("Scattering cross section is not a 2-D array.")

        if scattering.shape[0] != self._num_energy_groups or \
            scattering.shape[1] != self._num_energy_groups:
            raise ValueError("Size of scattering cross section invalid.")

        self._scattering = scattering

    @property
    def nu_fission(self):

        """The nu-fission macroscopic cross section.
        """

        return self._nu_fission

    @nu_fission.setter
    def nu_fission(self, nu_fission):

        """Sets the nu-fission macroscopic cross section.
        """

        nu_fission = np.array(nu_fission)

        if nu_fission.shape[0] != self._num_energy_groups:
            raise ValueError("Size of nu-fission cross section invalid.")

        self._nu_fission = nu_fission

    @property
    def diffusion(self):

        """The diffusion coefficient.
        """

        return self._diffusion

    @diffusion.setter
    def diffusion(self, diffusion):

        """Sets the diffusion coefficient.
        """

        diffusion = np.array(diffusion)

        if diffusion.shape[0] != self._num_energy_groups:
            raise ValueError("Size of diffusion coefficients invalid.")

        self._diffusion = diffusion

    @property
    def chi(self):

        """The fission energy spectrum.
        """

        return self._chi

    @chi.setter
    def chi(self, chi):

        """Sets the fission energy spectrum.
        """

        chi = np.array(chi)

        if chi.shape[0] != self._num_energy_groups:
            raise ValueError("Size of fission spectrum invalid.")

        self._chi = chi
