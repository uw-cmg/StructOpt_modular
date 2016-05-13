from ase.calculators.neighborlist import NeighborList

from structopt.tools import root, single_core, parallel


class HardSphereCutoff(object):
    """A relaxation module to ensure atoms in an individual are not too close together.
    This is often a preliminary relaxation before LAMMPS for VASP to ensure the models do not explode.
    """

    @single_core
    def __init__(self, cutoff=0.7):
        self.cutoff = cutoff


    @single_core
    def relax(individual):
        """Relaxes the individual using a hard-sphere cutoff method.
        Args:
            individual (Individual):  the individual to relax
        """
        radii = [2.0 for atom in individual]
        nl = NeighborList(radii, bothways=True, self_interaction=False)
        nl.update(individual)

        modified = True
        while modified:
            modified = False
            for atom in individual:
                indices, offsets = nl.get_neighbors(atom.index)
                for neigh in indices:
                    if individual.get_distance(atom.index, neigh) < self.cutoff:
                        individual.set_distance(atom.index, neigh, self.cutoff, fix=0.5)
                        modified = True
            nl.update(individual)

