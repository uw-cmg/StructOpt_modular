#from structopt.common.individual import Individual
from structopt.cluster import Cluster
from structopt.tools.dictionaryobject import DictionaryObject

import os
import logging
logging.parameters = DictionaryObject({
    "path": os.path.join(os.getcwd(), 'lammps_results'),
    "rank": 0,
    "ncores": 1
})

generator_params = DictionaryObject({
                       "read_xyz": {"filename": "ZrCu_xtal.xyz"}
                   })

relaxations = DictionaryObject({
    "LAMMPS": {"order": 0,
               "kwargs": {"use_mpi4py": True,
                      "MPMD": 0,
                      "keep_files": True,
                      "min_style": "cg",
                          "min_modify": "line quadratic",
#                      "minimize": "1e-8 1e-8 5000 10000",
                      "minimize": "1e-8 1e-8 1 1",
                      "pair_style": "eam/alloy",
                      "potential_file": "$STRUCTOPT_HOME/potentials/ZrCu.lammps.eam",
                      "thermo_steps": 0}}
})


individual = Cluster(id=None,
           load_modules=True,
           relaxation_parameters=relaxations, fitness_parameters=None,
           mutation_parameters=None,
           pso_moves_parameters=None,
           generator_parameters=generator_params)

individual.relax()

print(individual.LAMMPS)
