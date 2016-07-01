Parameters
##########


Global Parameters
=================

Seed

Initial Structures
==================

Crossovers
==========

Selections
==========

Mutations
=========

Relaxations
===========

The names of the modules to be used must be in a list at the top level of the relaxations parameters. For each module, a subsection at the top level of the relaxation paramters is necessary to specify the parameters for each module. The modules all need different paramters (refer tot individual modules sections below).

In addition to the module-specific parameters, each module requires two parallelization entries: ``use_mpi4py`` and ``MPMD_cores_per_structure``. These two entries are mutually exclusive, meaning that only one can be turned on at a time. ``use_mpi4py`` can take two values, ``true`` or ``false`` depending on whether the module should use the `one-structure-per-core <>`_ parallelization.

``MPMD_cores_per_structure`` should be set to ``0`` to disable it (i.e. if ``use_mpi4py`` is ``true``), but otherwise specifies the number of cores that each process/structure should be allocated within the ``MPI_Comm_spawn_multiple`` command. There are two types of valid values for this parameter: 1) an integer, or 2) a string of two integers specifying the minimum and maximum number of cores allowed (with syntax ``"4-16"``). ``MPMD_cores_per_structure`` can also take the value of ``"any"``, and StructOpt will use as many cores as it can to run each individual.

Example::

    "relaxations": {
        "modules": [
            "LAMMPS"
        ],
        "LAMMPS": {
            "use_mpi4py": true,
            "MPMD_cores_per_structure": 0,
            "keep_files": true,
            "min_style": "cg\nmin_modify line quadratic",
            "minimize": "1e-8 1e-8 5000 10000",
            "pair_style": "eam/alloy",
            "potential_file": "$STRUCTOPT_HOME/potentials/ZrCuAl2011.eam.alloy",
            "thermo_steps": 1000
        }
    }


LAMMPS
++++++

VASP
++++

Fitnesses
=========

See the Relaxations section above.

FEMSIM
++++++

STEM
++++

Predators
=========

Fingerprinters
==============

