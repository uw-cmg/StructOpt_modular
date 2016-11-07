import random
from ase import Atom, Atoms

from structopt.common.crossmodule import get_particle_radius

def ellipsoid(atomlist, fill_factor=0.74, radii=None, ratio=[1, 1, 1], cell=None):
    """Generates a random ellipsoid by rejection sampling. Not the most efficient
    code but good at garaunteeing an even distribution inside the ellipsoid.

    Parameters
    ----------
    fill_factor : float
        Determines how "close" the atoms are packed. 
        See structopt.tools.get_particle_radius for description
    radii : list (3 elements)
        The size, in angstroms, of the ellipsoid in the x, y and z
        direction
    ratio : list (3 elements)
        The ratio of the dimensions of the ellipsoid in the x, y
        and z direction.
    cell : list (3  elements)
        The size, in angstroms, of the dimensions that holds the
        atoms object
    """
    

    # Get the dimensions of the ellipsoid if required
    if radii is None:
        radius = get_particle_radius(atomlist, fill_factor)
        a = radius / ((ratio[1]/ratio[0])**(1.0/3.0) * (ratio[2]/ratio[0])**(1.0/3.0))
        b = radius / ((ratio[0]/ratio[1])**(1.0/3.0) * (ratio[2]/ratio[1])**(1.0/3.0))
        c = radius / ((ratio[0]/ratio[2])**(1.0/3.0) * (ratio[1]/ratio[2])**(1.0/3.0))
    else:
        a, b, c = radii

    # Create a list of random order of the atoms
    chemical_symbols = []
    for atom in atomlist:
        chemical_symbols += [atom[0]] * atom[1]
    random.shuffle(chemical_symbols)

    # Add atoms iteratively only if they fall inside the ellipsoid
    indiv = Atoms()
    while len(chemical_symbols) > 0:
        x, y, z = random.uniform(-a, a), random.uniform(-b, b), random.uniform(-c, c)
        if x**2/a**2 + y**2/b**2 + z**2/c**2 <= 1:
            indiv.extend(Atom(chemical_symbols.pop(), [x, y, z]))

    if cell is not None:
        indiv.set_cell(cell)
        indiv.center()
        indiv.set_pbc(True)

    return indiv
