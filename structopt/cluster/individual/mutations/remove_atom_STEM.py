import random

import numpy as np
from scipy.ndimage import filters

from structopt.tools import CoordinationNumbers
from structopt.tools import get_avg_radii
from structopt.common.individual.fitnesses import STEM

def remove_atom_STEM(individual, STEM_parameters,
                     filter_size=1, remove_CN=11, remove_cutoff=0.5, max_cutoff=0.5):

    """Moves surface atoms around based on the difference in the target
    and individual STEM image

    Parameters
    ----------
    STEM_parameters : dict
        Parameters for the STEM calculation. Ideally should be the same as the ones
        used for the STEM fitness/relaxation
    remove_CN : int
        The maximum coordination number considered as a surface atom to remove.    
    surf_CN : int
        The maximum coordination number considered as a surface atom to move an
        atom to
    filter_size : float
        Filter size for choosing local maximum in the picture. Filter size is equal
        to average_bond_length * resolution * filter_size.
    move_cutoff : float
        The search radius for selecting an atom to move near a high intensity point.
        Defaults to the average bond distance
    surf_cutoff : float
        The search radius for selecting a surface site near a low intensity point
        Defaults to the average bond distance
    """

    module = STEM(STEM_parameters)
    module.generate_target()
    target = module.target

    image, x_shift, y_shift = module.cross_correlate(module.get_image(individual))
    contrast = image - target
    max_max = np.max(contrast)

    # Determine filter size for locating local minimum
    cutoff = get_avg_radii(individual) * 2 * 1.1
    remove_cutoff *= cutoff
    resolution = module.parameters['resolution']
    size = cutoff * resolution * filter_size    

    ###################################
    ## Code for testing the contrast ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(contrast, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * 10))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * 10))
    # plt.show()
    # import sys; sys.exit()

    data_max = filters.maximum_filter(contrast, size=size)
    maxima = ((contrast == data_max) & (contrast > max_max * max_cutoff))
    if len(maxima) == 0:
        return False
    max_coords = np.argwhere(maxima)
    max_xys = (max_coords[:,::-1] - np.array([[x_shift, y_shift]])) / resolution
    max_intensities = np.asarray([data_max[tuple(coord)] for coord in max_coords])
    max_intensities /= sum(max_intensities)

    ###################################
    ## Code for testing the max find ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(maxima, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * 10))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * 10))
    # plt.show()
    # print(len(max_intensities))
    # import sys; sys.exit()

    # Get indices of atoms considered to be moved and sites to move too
    CNs = CoordinationNumbers(individual)
    positions = individual.get_positions()

    remove_indices = [i for i, CN in enumerate(CNs) if CN <= remove_CN]
    remove_xys = positions[list(remove_indices)][:, :2]

    # Randomly choose local maxima and minima locations from contrast weighted
    # by their intensity
    high_xy_index = np.random.choice(np.arange(len(max_xys)), p=max_intensities)
    high_xy = max_xys[high_xy_index]

    # Choose move_atom (surf_atom) from within the move_cutoff (surf_cutoff)
    # of high_xy (low_xy)
    dists_remove_xys = np.linalg.norm(remove_xys - high_xy, axis=1)
    indices_remove_xys = [i for i, d in zip(remove_indices, dists_remove_xys) if d < remove_cutoff]

    ########################
    ## Test atoms to move ##
    ########################
    # from ase.visualize import view
    # for i in indices_move_xys:
    #     individual[i].symbol = 'Mo'
    # view(individual)
    # import sys; sys.exit()

    if len(indices_remove_xys) == 0:
        remove_index = np.argmin(dists_remove_xys)
    else:
        remove_index = random.choice(indices_remove_xys)

    individual.pop(remove_index)

    return
