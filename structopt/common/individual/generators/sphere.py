import json
import random
import numpy as np
from ase import Atom, Atoms
from ase.visualize import view
from ase.data import atomic_numbers, reference_states

class Sphere(object):
    """Class that handles the generation of a sphere."""

    def __init__(self, atomlist, **kwargs):
        self.atomlist = atomlist
        self.fill_factor = 0.7
        self.radius = None
        self.min_dist_factor = 0.7
        self.tries_b4_expand = 100
        
        for kw in kwargs:
            setattr(self, kw, kwargs[kw])
    
    def random_three_vector(self):
        """Generates a random 3D unit vector (direction) with a 
        uniform spherical distribution
        Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
        :return:
        """

        phi = np.random.uniform(0,np.pi*2)
        costheta = np.random.uniform(-1,1)

        theta = np.arccos( costheta )
        x = np.sin( theta) * np.cos( phi )
        y = np.sin( theta) * np.sin( phi )
        z = np.cos( theta )
        return (x,y,z)

    def get_avg_radii(self):
        """Returns the average atomic radius of a list of
        atoms. The radius is the radius of the close packed sphere
        in a given crystal structure
        """

        atomlist = self.atomlist
    
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

    def get_particle_radius(self):
        """Returns an estimated nanoparticle radius given a
        concentration of atoms and void fraction of the sphere.
        Given an average sphere, this is given by the formula

        R_sphere = (n_tot / f)**(1.0/3.0) * R_atom

        where n_tot is the total number of atoms and f is
        the fill factor of the particle.
        """

        atomlist = self.atomlist
        fill_factor = self.fill_factor
        
        n_tot = sum([atom[1] for atom in atomlist])
        R_atom = self.get_avg_radii()
        R_sphere = (n_tot / fill_factor)**(1.0/3.0) * R_atom

        return R_sphere

    def generate(self):
        """Generates a random sphere of particles given an
        atomlist and radius. If radius is None, one is 
        automatically estimated. min_dist and tries_b4_expand
        are parameters that govern how stricly the proximity
        of atoms are enforced.
        """

        atomlist = self.atomlist
        radius = self.radius
        fill_factor = self.fill_factor
        min_dist_factor = self.min_dist_factor
        tries_b4_expand = self.tries_b4_expand

        indiv = Atoms()
        if hasattr(self, 'cell'):
            indiv.set_cell(self.cell)
            indiv.set_pbc(True)

        if radius is None:
            radius = self.get_particle_radius()

        avg_radii = self.get_avg_radii()
        min_dist = min_dist_factor * 2 * avg_radii

        # Create a list of random order of the atoms
        all_atoms = []
        for atom in atomlist:
            all_atoms += [atom[0]] * atom[1]

        random.shuffle(all_atoms)

        for i, atom in enumerate(all_atoms):
            # Always add the first atom to the origin
            if i == 0:
                indiv.append(Atom(atom, (0, 0, 0)))            
                continue

            # We want to assure atoms are not close together. We do this
            # by trial and error. At a certain point, we expand the radius
            # and try again. This is 
            j = 0
            while j <= tries_b4_expand:

                # Get a random coordinate in the sphere
                unit_vec = np.array(self.random_three_vector())
                D = radius * random.uniform(0, 1.0) ** (1.0/3.0)
                coord = unit_vec * D

                # Check the distances with other atoms
                pos = indiv.get_positions()
                dists = [np.linalg.norm(xyz - coord) for xyz in pos]
                if min(dists) < min_dist and j < tries_b4_expand:
                    j += 1
                elif j == tries_b4_expand:
                    print('expanded')
                    radius *= 1.1
                    j = 0
                else:
                    break

            indiv.append(Atom(atom, coord))

        return indiv
