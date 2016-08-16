import numpy as np
import scipy.ndimage.filters as filters
from scipy.optimize import fmin, brute

import structopt.common.individual.fitnesses
from structopt.tools.analysis import NeighborList
from structopt.tools import root, single_core, parallel
from structopt.tools import get_avg_radii, rotation_matrix

class STEM(structopt.common.individual.fitnesses.STEM):
    """Rotates the individual to obtain a better match with STEM image.
    This is done by taking one bright column in the bulk region and reading
    its nearest neighbors as a projection into the xy plane. A bulk atom close
    to the individual COM and its nearest neighbors in the individual is selected
    and rotated to optimize the matching of its projection in the xy plane
    with that of the STEM image.

    Parameters
    ---------
    HWHM : float
        The half-width half-maximum of the gaussian function used for the 
        point spread function.
    dimensions : list
        The x and y dimensions of STEM image.
    resolution : float
        The pixels per angstrom resolution
    gridsize : int
        Given a individual nearest neighbor unit to be optimized with the
        STEM image, the gridsize determines how fine the search for
        the rotation be. Tests indicate a gridsize of 10 is suitable.
    """

    def __init__(self, parameters={}):
        if 'rotation_grid' not in parameters:
            parameters.setdefault('rotation_grid', 10)
        super().__init__(parameters)

    def relax(self, individual):
        if self.psf is None:
            self.generate_psf()
        if self.target is None:
            self.generate_target()

        steps = self.parameters['rotation_grid']
        bonds = self.get_bulk_bonds(individual)
        projection = self.get_STEM_projection(individual)
        solution = brute(self.epsilon,
                         args=(bonds, projection),
                         ranges=(slice(0, np.pi, np.pi/steps),
                                 slice(-1, 1, 2/steps),
                                 slice(0, 2*np.pi, 2*np.pi/steps)),
                         finish=fmin)

        phi, costheta, a = solution
        theta = np.arccos(costheta)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        individual.rotate([x, y, z], a, center='COP')
        individual.STEM = self.fitness(individual)

        return

    @staticmethod
    def epsilon(rotation, bonds, projection):
        """Calculates the difference in the projected xy coordinates
        of a rotated bonding center of that of the individual
        and of the STEM image

        Parameters
        ----------
        rotation : list [phi, cos(theta), a]
            Contains phi (rotation around z-axis) and costheta (rotation 
            along the z-axis) and the angle to rotate the atoms object
        bonds : list of [vx, vy, vz] vectors
            The bonds around an atomic center in the bulk of the individual.
        projection : list of [vx, vy] vectors
            The bonds around an atomic center from the STEM image

        Output
        ------
        out : float
            The sum squared difference between projected [vx, vy] vectors
            after performing rotation.
        """

        phi, costheta, a = rotation
        theta = np.arccos(costheta)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        rotate = rotation_matrix(np.array([x, y, z]), a)
        bonds = np.asarray([[np.dot(rotate, v)[:2] for v in bonds]])
        projection = np.asarray([projection])
        dists = np.linalg.norm(np.transpose(bonds, (1, 0, 2)) - projection, axis=2)
        total_error = np.sum(np.square(np.amin(dists, axis=1)))

        return total_error

    def get_bulk_bonds(self, individual):
        """Returns the bonds to a bulk atom near the center of mass"""

        # Get a bulk atom near the center of the particle
        pos = individual.get_positions()
        com = individual.get_center_of_mass()
        dists_from_com = np.linalg.norm(pos - com, axis=1)
        bulk_atom_index = np.argmin(dists_from_com)
        bulk_atom_pos = pos[bulk_atom_index]

        # Get the neighbors of the bulk atom
        NNs = NeighborList(individual)
        bonds = np.asarray([pos[i] for i in NNs[bulk_atom_index]]) - bulk_atom_pos

        return bonds

    def get_STEM_projection(self, individual):
        """Gets the projection of bonds from the brighest STEM image"""

        if self.target is None:
            self.generate_target()

        parameters = self.parameters
        target = self.target

        # Get a list of xy positions from analyzing local maxima in STEM image
        # as well as the position of the brighest atom
        chemical_symbols = individual.get_chemical_symbols()
        unique_symbols = set(chemical_symbols)
        atomlist = [[symbol, chemical_symbols.count(symbol)]
                    for symbol in unique_symbols]
        cutoff = get_avg_radii(atomlist) * 2 * 1.1
        size = cutoff / 8 * parameters['resolution']

        data_max = filters.maximum_filter(target, size=size)
        maxima = ((target == data_max) & (target > 0.1)) # Filter out low maxima
        max_pixel = np.argmax(target)
        y, x = np.unravel_index(np.argmax(target), target.shape)
        max_pos = np.array([x, y]) / np.asarray(parameters['resolution'])
        pos = np.asarray(np.argwhere(maxima)[:,::-1] / parameters['resolution'])
        vecs = pos - max_pos
        dists = np.linalg.norm(vecs, axis=1)
        vecs = vecs[((dists < cutoff) & (dists > 0))]

        return vecs
