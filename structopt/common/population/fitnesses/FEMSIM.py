import logging
import subprocess
import math
import time
from collections import defaultdict
from mpi4py import MPI

import structopt
from structopt.tools import root, single_core, parallel


@root
def fitness(population):
    """Perform the FEMSIM fitness calculation on an entire population.

    Args:
        population (Population): the population to evaluate
    """
    to_fit = [individual for individual in population if not individual._fitted]

    if to_fit:
        cores_per_individual = structopt.parameters.globals.ncores // len(to_fit)
        # Round cores_per_individual down to nearest power of 2
        if cores_per_individual == 0:
            # We have more individuals than cores, so each fitness scheme needs to be run multiple times
            cores_per_individual = 1

        pow(2.0, math.floor(math.log2(cores_per_individual)))

        logger = logging.getLogger('output')

        # Setup each individual and get the inputs for each individual that need to be passed into the spawn
        multiple_spawn_args = defaultdict(list)
        for individual in to_fit:
            spawn_args = individual.fitnesses.FEMSIM.get_spawn_args(individual)
            for arg_name, arg in spawn_args.items():
                multiple_spawn_args[arg_name].append(arg)

        # Make sure each key in `multiple_spawn_args` has the same number of elements and set `count` to that value
        count = -1
        for key, value in multiple_spawn_args.items():
            if count != -1:
                assert len(value) == count
            count = len(value)

        # Create MPI.Info objects from the kwargs dicts in multiple_spawn_args['info'] for each rank
        infos = [MPI.Info.Create() for _ in multiple_spawn_args['info']]
        for i, info in enumerate(infos):
            for key, value in multiple_spawn_args['info'][i].items():
                info.Set(key, value)

        # Run the multiple spawn
        print("Spawning {} femsim processes, each with {} cores".format(count, cores_per_individual))
        MPI.COMM_SELF.Spawn_multiple(command=multiple_spawn_args['command'],
                                     args=multiple_spawn_args['args'],
                                     maxprocs=[cores_per_individual]*count,
                                     info=infos
                                     )

        # Collect the results for each chisq and return them
        for i, individual in enumerate(population):
            while not individual.fitnesses.FEMSIM.has_finished():
                time.sleep(0.1)
            time.sleep(0.1)
            chisq = individual.fitnesses.FEMSIM.get_chisq(individual)
            individual.FEMSIM = chisq
            logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, chisq))

    return [individual.FEMSIM for individual in population]

