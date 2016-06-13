Setup
#####


Setup an Anaconda virtual environment with the latest version of Python 3.


mpi4py
======

On Madion's ACI:

::

    module load compile/intel
    module load mpi/intel/openmpi-1.10.2

Follow these instructions:  http://mpi4py.readthedocs.io/en/stable/install.html#using-distutils

::

    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-X.Y.tar.gz
    tar -zxf mpi4py-X.Y.tar.gz
    cd mpi4py-X.Y
    python setup.py build
    python setup.py install --user

Try at least the first test:  http://mpi4py.readthedocs.io/en/stable/install.html#testing


Python Libraries
================

::

    pip install ase
    pip install natsorted
