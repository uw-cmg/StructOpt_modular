Parallelism
###########

Note: mpi4py needs to be installed from source. I followed the instructions here https://media.readthedocs.org/pdf/mpi4py/latest/mpi4py.pdf under "3.3: Using distutils".

::

    # Setup modules so that `mpi/intel/openmpi` is loaded and `mpiexec` finds that executable
    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-X.Y.tar.gz
    tar -zxf mpi4py-X.Y.tar.gz
    cd mpi4py-X.Y
    python setup.py build
    python setup.py install --user


mpi4py
======

Main idea:  One structure per core, or multiple structures per core that execute serially in a for-loop.

`mpi4py <https://mpi4py.scipy.org/docs/usrman/tutorial.html>`_ allows MPI commands to be run within python. Some fitness and relaxations modules may want to parallelize along the lines of one structure per core. mpi4py is used for this functionality.

Cases
-----

* ``ncores == len(population)``

* ``ncores < len(population)``

* ``ncores > len(population)`` (unused cores that were, for example, allocated for use within MDMP)

MPMD
====

Main idea:  Multiple cores per structure.

Multiple program, multiple data (MPMD) is a form of MPI parallelization where multiple commands can be passed to ``mpiexec`` separated by ``:``s. Each command is preceded by the ``-n`` option whcih specifies the number of cores to be used for that executable. The executable will then need to implement MPMD (often by splitting the world communicator up by color) (see `here <https://github.com/jjmaldonis/mpi-parallelization/blob/master/testmpi.f90>`_ for an example).
