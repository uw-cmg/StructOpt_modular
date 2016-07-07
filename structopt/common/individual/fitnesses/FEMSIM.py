import os
import logging
import numpy as np
import shutil

import structopt
from structopt.io import write_xyz
from structopt.tools import root, single_core, parallel


class FEMSIM(object):
    """Contains parameters and functions for running FEMSIM through Python."""

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.k = None
        self.parameters = self.read_inputs(parameters)
        self.vk = np.multiply(self.parameters.thickness_scaling_factor, self.vk)  # Multiply the experimental data by the thickness scaling factor

        # These parameteres do not need to exist between generations
        # They are used for before/after femsim processing
        self.base = None
        self.folder = None
        self.paramfilename = None


    @single_core
    def read_inputs(self, parameters):
        data = open(parameters.vk_data_filename).readlines()
        data.pop(0)  # Comment line
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        k, vk = zip(*data)
        # Set k and vk data for chi2 comparison
        self.k = np.array(k)
        self.vk = np.array(vk)
        return parameters


    @single_core
    def update_parameters(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self.parameters, key, value)


    @single_core
    def get_spawn_args(self, individual):
        """Returns a dictionary of arguments to be passed to MPI.COMM_SELF.Spawn which will be collected for all
        structures and concatenated into MPI.COMM_SELF.Spawn_multiple: 
        https://github.com/mpi4py/mpi4py/blob/2acfc552c42846628304e54a3b87e2bf3a59af07/src/mpi4py/MPI/Comm.pyx#L1555
        """
        self.setup_individual_evaluation(individual)
        femsim_command = os.environ['FEMSIM_COMMAND']
        args = [self.base, self.paramfilename]
        info = {'wdir': self.folder}
        return {'command': femsim_command, 'args': args, 'info': info}


    @single_core
    def get_chisq(self, individual):
        vk = self.get_vk_data()
        return self.chi2(vk)


    @single_core
    def setup_individual_evaluation(self, individual):

        logger = logging.getLogger('by-rank')

        logger.info('Received individual HI = {0} for FEMSIM evaluation'.format(individual.index))

        # Make individual folder and copy files there
        self.folder = os.path.abspath('Output-rank0/FEMSIMFiles/Individual{i}'.format(i=individual.index))
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)
        if not os.path.isfile(os.path.join(self.folder, self.parameters.vk_data_filename)):
            shutil.copy(self.parameters.vk_data_filename, os.path.join(self.folder, self.parameters.vk_data_filename))

        self.paramfilename = os.path.join(self.folder, "femsim.{}.in".format(individual.index))
        shutil.copy(self.parameters.parameter_filename, self.paramfilename)
        self.write_paramfile(individual)

        base = 'indiv{i}'.format(i=individual.index) # TODO Add generation number
        self.base = base


    @single_core
    def write_paramfile(self, individual):
        # Write structure file to disk so that the fortran femsim can read it in
        #ase.io.write('structure_{i}.xyz'.format(i=individual.index), individual)
        comment = "{} {} {}".format(self.parameters.xsize, self.parameters.ysize, self.parameters.zsize)
        write_xyz('structure_{i}.xyz'.format(i=individual.index), individual, comment=comment)

        with open(self.paramfilename, 'w') as f:
            f.write('# Parameter file for generation {gen}, individual {i}\n'.format(gen=None, i=individual.index))  # TODO Add generation number
            f.write('{}\n'.format(os.path.join(os.getcwd(), 'structure_{i}.xyz'.format(i=individual.index))))
            f.write('{}\n'.format(self.parameters.vk_data_filename))
            f.write('{}\n'.format(self.parameters.Q))
            f.write('{} {} {}\n'.format(self.parameters.nphi, self.parameters.npsi, self.parameters.ntheta))


    @single_core
    def has_finished(self):
        return os.path.exists(os.path.join(self.folder, 'vk_initial_{base}.txt'.format(base=self.base)))


    @single_core
    def get_vk_data(self):
        filename = os.path.join(self.folder, 'vk_initial_{base}.txt'.format(base=self.base))
        data = open(filename).readlines()
        os.remove(filename)
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        vk = np.array([vk for k, vk in data])
        return vk


    @single_core
    def chi2(self, vk):
        return np.sum(((self.vk - vk) / self.vk)**2) / len(self.k)

