import json
import random
import numpy as np
from ase import Atom, Atoms
from ase.visualize import view
from ase.data import atomic_numbers, reference_states

def random_three_vector():
    """Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """

    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return (x,y,z)

def get_avg_radii(atomlist):
    """Returns the average atomic radius of a list of
    atoms. The radius is the radius of the close packed sphere
    in a given crystal structure
    """

    # Get the average atomic radii of close packed atoms
    n_tot = sum([atom[1] for atom in atomlist])
    r = 0
    for atom in atomlist:
        n = atom[1]
        conc = float(n)/float(n_tot)
        atomic_number = atomic_numbers[atom[0]]
        struct = reference_states[atomic_number]['symmetry']
        if struct == 'fcc':
            a = reference_states[atomic_number]['a']
            r += conc * np.linalg.norm([a, a]) / 4.0
        elif struct == 'bcc':
            a = reference_states[atomic_number]['a']
            r += conc * np.linalg.norm([a, a, a]) / 4.0
        else:
            raise IOError('{} structure not supported yet'.format(struct))

    return r

def get_particle_radius(atomlist, fill_factor):
    """Returns an estimated nanoparticle radius given a
    concentration of atoms and void fraction of the sphere.
    Given an average sphere, this is given by the formula

    R_sphere = (n_tot / f)**(1.0/3.0) * R_atom

    where n_tot is the total number of atoms and f is
    the fill factor of the particle.
    """

    n_tot = sum([atom[1] for atom in atomlist])
    R_atom = get_avg_radii(atomlist)
    R_sphere = (n_tot / fill_factor)**(1.0/3.0) * R_atom

    return R_sphere

def sphere(atomlist, fill_factor=0.7, radius=None, cell=None):
    """Generates a random sphere of particles given an
    atomlist and radius. If radius is None, one is 
    automatically estimated. min_dist and tries_b4_expand
    are parameters that govern how stricly the proximity
    of atoms are enforced.
    """

    if not radius:
        radius = get_particle_radius(atomlist, fill_factor)

    avg_radii = get_avg_radii(atomlist)

    # Create a list of random order of the atoms
    chemical_symbols = []
    for atom in atomlist:
        chemical_symbols += [atom[0]] * atom[1]

    random.shuffle(chemical_symbols)

    unit_vec = np.array([np.array(random_three_vector()) for i in range(len(chemical_symbols))])
    D = radius * np.random.sample(size=len(chemical_symbols)) ** (1.0/3.0)
    positions = np.array([D]).T * unit_vec

    indiv = Atoms(symbols=chemical_symbols, positions=positions)

    if cell is not None:
        indiv.set_cell(cell)
        cell_center = np.sum(indiv.get_cell(), axis=1) / 2.0
        indiv.translate(indiv.get_center_of_mass() + cell_center)
        indiv.set_pbc(True)

    return indiv
