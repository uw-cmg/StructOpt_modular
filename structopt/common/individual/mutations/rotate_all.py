import random
import numpy as np

from structopt.tools import random_three_vector


def rotate_all(atoms, vector=None, angle=None, center=None):
    """Rotate all atoms around a single point. Most suitable for
    cluster calculations.

    Parameters
    ---------
    individual : StructOpt individual object or ase atoms
        StructOpt Individual or ase Atoms object to be rotated.
    angle : string or list
        A list of angles that will be chosen to rotate. If None,
        is randomly generated. Angle must be given in radians.
        If 'random' in list, a random angle is included.
    vector : string or list
        The list of axes in which to rotate the atoms around. If 
        None, is a randomly chosen direction. If 'random' in list,
        a random vector can be chosen.
    center : str or xyz iterable
        The center in which to rotate the atoms around. If None,
        defaults to center of mass. Acceptable strings are
        COM = center of mass
        COP = center of positions
        COU = center of cell

    Returns
    -------
    out: None
        Modifies the individual in place
    """

    # Initialize variables for ase.Atoms.rotate
    if angle is None:
        angle = random.uniform(30, 180) * np.pi / 180.0
    elif hasattr(angle, '__iter__'):
        angle = random.choice(angle)
        if angle == 'random':
            angle = random.uniform(30, 180) * np.pi / 180
        
    if vector is None:
        vector = random_three_vector()
    elif hasattr(vector, '__iter__'):
        vector = random.choice(vector)
        if angle == 'random':
            vector = random_three_vector()

    if center is None:
        center = 'COM'

    # Perform the rotation
    atoms.rotate(v=vector, a=angle, center=center)

    return vector, angle
