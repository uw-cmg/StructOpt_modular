import random

import numpy as np
from scipy.ndimage import filters

from structopt.tools import CoordinationNumbers
from structopt.tools import get_avg_radii
from structopt.common.individual.fitnesses import STEM

from ase.io import write

def move_surface_STEM(individual, STEM_parameters, surf_CN=10, filter_size=1):
    """Moves surface atoms around based on the difference in the target
    and individual STEM image

    Parameters
    ----------
    STEM_parameters : dict
        Parameters for the STEM calculation. Ideally should be the same as the ones
        used for the STEM fitness/relaxation
    surf_CN : int
        The maximum coordination number considered as a surface atom. Surface atoms
        are considered for choosing which atoms to move and where to move them
    filter_size : float
        Filter size for choosing local maximum in the picture. Filter size is equal
        to average_bond_length * resolution * filter_size.
    """

    module = STEM(STEM_parameters)
    module.generate_target()
    target = module.target

    image, x_shift, y_shift = module.cross_correlate(module.get_image(individual))
    contrast = image - target

    # Find a list of local maximum and local minimum in the image
    cutoff = get_avg_radii(individual) * 2 * 1.1
    resolution = module.parameters['resolution']        
    size = cutoff * resolution * filter_size

    data_max = filters.maximum_filter(contrast, size=size)
    maxima = ((contrast == data_max) & (contrast > 0.1)) # Filter out low maxima
    max_coords = np.argwhere(maxima)
    max_xys = (max_coords[:,::-1] + [x_shift, y_shift]) / resolution
    max_intensities = np.array([data_max[tuple(coord)] for coord in max_coords])

    data_min = filters.minimum_filter(contrast, size=size)
    minima = ((contrast == data_min) & (contrast < -0.1)) # Filter out high minim
    min_coords = np.argwhere(minima)
    min_xys = (min_coords[:,::-1] + [x_shift, y_shift]) / resolution
    min_intensities = np.absolute(np.array([data_min[tuple(coord)] for coord in min_coords]))

    # Get a list of surface atoms
    CNs = CoordinationNumbers(individual)
    positions = individual.get_positions()

    surf_CN = 11
    surf_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN]
    surf_positions = positions[list(surf_indices)]
    surf_xys = surf_positions[:, :2]

    # Randomly choose local maxima and minima locations from contrast
    high_xy_index = np.random.choice(np.arange(len(max_xys)),
                                     p=max_intensities/sum(max_intensities))
    low_xy_index = np.random.choice(np.arange(len(min_xys)),
                                    p=min_intensities/sum(min_intensities))
    high_xy = max_xys[high_xy_index]
    low_xy = min_xys[low_xy_index]

    # Choose atom nearest to the high_xy and position nearest to the low_xy
    move_index = np.argmin(np.linalg.norm(surf_xys - high_xy, axis=1))
    new_position = surf_positions[np.argmin(np.linalg.norm(surf_xys - low_xy, axis=1))]
    new_position_vec = (new_position - individual.get_center_of_mass()) / np.linalg.norm(new_position - individual.get_center_of_mass())    
    new_position += new_position_vec * cutoff/2

    positions[surf_indices[move_index]] = new_position
    individual.set_positions(positions)

    return
