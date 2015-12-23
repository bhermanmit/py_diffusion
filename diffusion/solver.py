"""Solver module.
"""

from scipy.sparse import coo_matrix
import scipy.sparse.linalg
import numpy as np


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

        self._num_mesh = self._mesh.part_mesh.sum()
        self._num_parts = self._mesh.num_parts
        self._num_energy_groups = self._materials.num_energy_groups

        self._mesh.generate()
        self._materials.generate()

        self._left_albedo = 0.0
        self._right_albedo = 0.0

        self._num_eigs = 1

    @property
    def num_eigs(self):

        """The number of eigenvalues to solve for.
        """

        return self._num_eigs

    @num_eigs.setter
    def num_eigs(self, num_eigs):

        """Sets the number eigenvalues
        """

        self._num_eigs = num_eigs

    def build_loss_matrix(self):

        """Builds left hand side loss matrix.
        """

        # compute number of non zero elements
        num_elements = ((self._num_mesh - 2)*3 + 4) * \
            self._num_energy_groups
        num_elements += self._num_energy_groups**2 * self._num_mesh - \
            self._num_energy_groups * self._num_mesh

        # allocate vectors
        rows = np.zeros((num_elements,), dtype=np.int)
        cols = np.zeros((num_elements,), dtype=np.int)
        vals = np.zeros((num_elements,))
        count = 0

        # loop around rows
        for i_row in xrange(self._num_mesh * self._num_energy_groups):

            # get indices
            mesh_index, group_index = self._matrix_to_indices(i_row)

            # look up material index
            material_index = self._mesh.get_material(mesh_index)

            # get material handle
            material = self._materials[material_index]

            # get data
            width_self = self._mesh.get_width(mesh_index)
            diffusion_self = material.diffusion[group_index]
            removal_self = material.removal[group_index]

            # leakage left side neighbor
            if mesh_index != 0:

                left_index = mesh_index - 1
                width_left = self._mesh.get_width(left_index)
                left_material_index = self._mesh.get_material(left_index)
                left_material = self._materials[material_index]

                diffusion_left = left_material.diffusion[group_index]

                val = -2.0*diffusion_left*diffusion_self / \
                    (diffusion_left*width_self + diffusion_self*width_left)

                # record in sparse matrix
                i_col = self._indices_to_matrix(left_index, group_index)
                rows[count] = i_row
                cols[count] = i_col
                vals[count] = val / width_self
                count += 1

            # leakage right side neighbor
            if mesh_index < self._num_mesh - 1:

                right_index = mesh_index + 1
                width_right = self._mesh.get_width(right_index)
                right_material_index = self._mesh.get_material(right_index)
                right_material = self._materials[material_index]

                diffusion_right = right_material.diffusion[group_index]

                val = -2.0*diffusion_right*diffusion_self / \
                    (diffusion_right*width_self + diffusion_self*width_right)

                # record in sparse matrix
                i_col = self._indices_to_matrix(right_index, group_index)
                rows[count] = i_row
                cols[count] = i_col
                vals[count] = val / width_self
                count += 1

            # handle self
            val = 0.0

            # leakage left
            if mesh_index == 0:

                val += 2.0*diffusion_self*(1.0 - self._left_albedo) / \
                    (4*diffusion_self*(1.0 + self._left_albedo) +
                     (1.0 - self._left_albedo)*width_self)

            else:

                val += 2.0*diffusion_left*diffusion_self / \
                    (diffusion_left*width_self + diffusion_self*width_left)

            # leakage right
            if mesh_index == self._num_mesh - 1:

                val += 2.0*diffusion_self*(1.0 - self._right_albedo) / \
                    (4*diffusion_self*(1.0 + self._right_albedo) +
                     (1.0 - self._right_albedo)*width_self)

            else:

                val += 2.0*diffusion_right*diffusion_self / \
                    (diffusion_right*width_self + diffusion_self*width_right)

            # divide leakage by width
            val /= width_self

            # remove cross section
            val += removal_self

            # record in sparse matrix
            rows[count] = i_row
            cols[count] = i_row
            vals[count] = val
            count += 1

            # scattering
            for h in xrange(self._num_energy_groups):

                # dont record within group scattering, removal has it
                if h == group_index:
                    continue

                # in-scattering cross section
                val = -material.scattering[h, group_index]

                # record in sparse matrix
                i_col = self._indices_to_matrix(mesh_index, h)
                rows[count] = i_row
                cols[count] = i_col
                vals[count] = val
                count += 1

        # construct sparse matrix
        dimension = self._num_mesh * self._num_energy_groups
        loss_matrix = coo_matrix((vals, (rows, cols)),
                                 shape=(dimension, dimension)).tocsc() 

        return loss_matrix

    def build_production_matrix(self):

        """Builds right hand side production matrix.
        """

        # compute number of non zero elements
        num_elements = self._num_mesh * self._num_energy_groups**2

        # allocate vectors
        rows = np.zeros((num_elements,), dtype=np.int)
        cols = np.zeros((num_elements,), dtype=np.int)
        vals = np.zeros((num_elements,))
        count = 0

        # loop around rows
        for i_row in xrange(self._num_mesh * self._num_energy_groups):

            # get indices
            mesh_index, group_index = self._matrix_to_indices(i_row)

            # look up material index
            material_index = self._mesh.get_material(mesh_index)

            # get material handle
            material = self._materials[material_index]

            # get data
            chi = material.chi[group_index]

            # nu fission
            for h in xrange(self._num_energy_groups):

                # nu fission cross section
                val = chi*material.nu_fission[h]

                # record in sparse matrix
                i_col = self._indices_to_matrix(mesh_index, h)
                rows[count] = i_row
                cols[count] = i_col
                vals[count] = val
                count += 1

        # construct sparse matrix
        dimension = self._num_mesh * self._num_energy_groups
        prod_matrix = coo_matrix((vals, (rows, cols)),
                                 shape=(dimension, dimension)).tocsc()

        return prod_matrix

    def solve(self):

        """Solves the diffusion system.
        """

        loss_matrix = self.build_loss_matrix()
        prod_matrix = self.build_production_matrix()

        eigs, vectors = scipy.sparse.linalg.eigs(A=prod_matrix, M=loss_matrix,
                                                 which='LR', k=self._num_eigs)

        eigs = eigs.real
        vectors = vectors.reshape((self._num_energy_groups, self._num_mesh,
                                   self._num_eigs))
        vectors = vectors.real

        self._eigs = eigs
        self._vectors = vectors

    def extract_eigenvalues(self):

        """Extracts the eigenvalues.
        """

        return self._eigs

    def extract_eigenvectors(self):

        """Extracts the eigenvectors.
        """

        slab_pos = self._mesh.extract_x()

        return slab_pos, self._vectors

    def _matrix_to_indices(self, matrix):

        """Converts matrix row/col to indices.
        """

        mesh_index = matrix % self._num_mesh
        group_index = (matrix % (self._num_mesh*self._num_energy_groups)) // \
            self._num_mesh

        return mesh_index, group_index

    def _indices_to_matrix(self, mesh_index, group_index):

        """Converst indices to matrix row/col.
        """

        return mesh_index + group_index*self._num_mesh
