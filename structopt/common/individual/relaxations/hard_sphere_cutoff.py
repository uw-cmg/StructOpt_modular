from ase.calculators.neighborlist import NeighborList

from structopt.tools import root, single_core, parallel


#class HardSphereCutoff(object):
class hard_sphere_cutoff(object):
    """A relaxation module to ensure atoms in an individual are not too close together.
    This is often a preliminary relaxation before LAMMPS for VASP to ensure the models do not explode.
    """

    @single_core
    def __init__(self, cutoff=0.7):
        self.cutoff = cutoff


    @single_core
    def relax(self, individual):
        """Relaxes the individual using a hard-sphere cutoff method.
        Args:
            individual (Individual):  the individual to relax
        """
        print("Relaxing individual {} with hard-sphere cutoff method".format(individual.index))
        radii = [2.0 for atom in individual]
        nl = NeighborList(radii, bothways=True, self_interaction=False)
        nl.update(individual)

        ntries = 0
        modified = True
        while modified and ntries < 100:
            modified = False
            for atom in individual:
                indices, offsets = nl.get_neighbors(atom.index)
                for neigh in indices:
                    if individual.get_distance(atom.index, neigh) < self.cutoff:
                        individual.set_distance(atom.index, neigh, self.cutoff, fix=0.5)
                        modified = True
            nl.update(individual)
            individual.wrap()
            ntries += 1
        if ntries == 100:
            print("WARNING! Iterated through the hard-sphere cutoff relaxation 100 times and it still did not converge!")

