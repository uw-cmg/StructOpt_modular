import random

import numpy as np
from scipy.ndimage import filters
from ase import Atom

from structopt.tools import CoordinationNumbers
from structopt.tools import get_avg_radii
from structopt.common.individual.fitnesses import STEM

def add_atom_STEM(individual, STEM_parameters, elements=None, p=None,
                  filter_size=1, surf_CN=11, surf_cutoff=0.5, min_cutoff=0.5):

    """Moves surface atoms around based on the difference in the target
    and individual STEM image

    Parameters
    ----------
    STEM_parameters : dict
        Parameters for the STEM calculation. Ideally should be the same as the ones
        used for the STEM fitness/relaxation
    move_CN : int
        The maximum coordination number considered as a surface atom to move.    
    surf_CN : int
        The maximum coordination number considered as a surface atom to move an
        atom to
    filter_size : float
        Filter size for choosing local maximum in the picture. Filter size is equal
        to average_bond_length * resolution * filter_size.
    surf_cutoff : float
        The search radius for selecting a surface site near a low intensity point
        Defaults to the average bond distance
    """

    module = STEM(STEM_parameters)
    module.generate_target()
    target = module.target

    image, x_shift, y_shift = module.cross_correlate(module.get_image(individual))
    contrast = image - target
    min_min = np.min(contrast)

    # Determine filter size for locating local minimum
    cutoff = get_avg_radii(individual) * 2 * 1.1
    surf_cutoff *= cutoff
    resolution = module.parameters['resolution']
    size = cutoff * resolution * filter_size    

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

    data_min = filters.minimum_filter(contrast, size=size)
    minima = ((contrast == data_min) & (contrast < min_min * min_cutoff))
    if len(minima) == 0:
        return False
    min_coords = np.argwhere(minima)
    min_xys = (min_coords[:,::-1] - [x_shift, y_shift]) / resolution
    min_intensities = np.asarray([data_min[tuple(coord)] for coord in min_coords])
    min_intensities = np.absolute(min_intensities)
    min_intensities /= sum(min_intensities)

    ###################################
    ## Code for testing the min find ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(data_min, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * 10))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * 10))
    # plt.show()
    # print(len(min_intensities))
    # import sys; sys.exit()

    # Get indices of atoms considered to be moved and sites to move too
    CNs = CoordinationNumbers(individual)
    positions = individual.get_positions()

    surf_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN]
    surf_positions = positions[list(surf_indices)]
    COM = surf_positions.mean(axis=0)
    vec = surf_positions - COM
    vec /= np.array([np.linalg.norm(vec, axis=1)]).T
    epsilon = np.array([np.random.random(len(surf_positions)) * 0.5 + 0.5]).T
    surf_positions = surf_positions + vec * cutoff * epsilon
    surf_xys = surf_positions[:, :2]

    # Randomly choose local maxima and minima locations from contrast weighted
    # by their intensity
    low_xy_index = np.random.choice(np.arange(len(min_xys)), p=min_intensities)
    low_xy = min_xys[low_xy_index]

    ########################
    ## Test atoms to move ##
    ########################
    # from ase.visualize import view
    # for i in indices_move_xys:
    #     individual[i].symbol = 'Mo'
    # view(individual)
    # import sys; sys.exit()

    dists_surf_xys = np.linalg.norm(surf_xys - low_xy, axis=1)
    indices_surf_xys = [i for i, d in enumerate(dists_surf_xys) if d < surf_cutoff]

    ########################
    ## Test atoms to move ##
    ########################
    # from ase import Atom
    # from ase.visualize import view
    # for i in indices_surf_xys:
    #     individual.append(Atom('Mo', surf_positions[i]))
    # view(individual)
    # import sys; sys.exit()

    # Pick the surface site to add the atom to
    if len(indices_surf_xys) == 0:
        surf_index = np.argmin(dists_surf_xys)
    else:
        surf_index = random.choice(indices_surf_xys)
    new_position = surf_positions[surf_index]

    # Choose the element to add
    if elements is None and p is None:
        syms = individual.get_chemical_symbols()
        elements = np.unique(syms)
        n = len(syms)
        p = [syms.count(element) / n for element in elements]
    element = np.random.choice(elements, p=p)

    individual.append(Atom(element, new_position))

    return
