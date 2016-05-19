Parallelism
###########

mpi4py
======

Main idea:  One structure per core, or multiple structures per core that execute serially in a for-loop.

`mpi4py <https://mpi4py.scipy.org/docs/usrman/tutorial.html>`_ allows MPI commands to be run within python. Some fitness and relaxations modules may want to parallelize along the lines of one structure per core. mpi4py is used for this functionality.

Cases
-----

1) ``ncores == len(population)``

2) ``ncores < len(population)``

3) ``ncores > len(population)`` (unused cores)

MPMD
====

Main idea:  Multiple cores per structure.

Multiple program, multiple data (MPMD) is a form of MPI parallelization where multiple commands can be passed to ``mpiexec`` separated by ``:``s. Each command is preceded by the ``-n`` option whcih specifies the number of cores to be used for that executable. The executable will then need to implement MPMD (often by splitting the world communicator up by color) (see `here <https://github.com/jjmaldonis/mpi-parallelization/blob/master/testmpi.f90>`_ for an example).
