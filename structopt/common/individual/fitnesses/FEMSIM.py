import os
import logging
import numpy as np
import shutil

import structopt
from structopt.io import write_xyz
from structopt.tools import root, single_core, parallel


class FEMSIM(object):
    """ """

    @single_core
    def __init__(self):
        self.parameters = self.read_inputs()

        self.vk = np.multiply(self.parameters.thickness_scaling_factor, self.vk)  # Multiply the experimental data by the thickness scaling factor


    @single_core
    def read_inputs(self):
        parameters = structopt.parameters.fitnesses.FEMSIM
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
    def get_command(self, individual):
        self.setup_individual_evaluation(individual)
        femsim_command = os.environ['FEMSIM_COMMAND']
        return '-wdir {folder} {femsim_command} {base} {paramfile}'.format(folder=self.folder, femsim_command=femsim_command, base=self.base, paramfile=self.paramfilename)


    @single_core
    def get_chisq(self, individual):
        vk = self.get_vk_data()
        return self.chi2(vk)


    @single_core
    def setup_individual_evaluation(self, individual):

        logger = logging.getLogger('by-rank')

        logger.info('Received individual HI = {0} for FEMSIM evaluation'.format(individual.index))

        # Make individual folder and copy files there
        self.folder = os.path.abspath('{filename}-rank0/FEMSIMFiles/Individual{i}'.format(filename=structopt.parameters.globals.output_filename, i=individual.index))
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
            f.write('{}\n'.format(self.parameters.thickness_scaling_factor))


    @single_core
    def get_vk_data(self):
        data = open(os.path.join(self.folder, 'vk_initial_{base}.txt'.format(base=self.base))).readlines()
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        vk = np.array([vk for k, vk in data])
        return vk


    @single_core
    def chi2(self, vk):
        return np.sum(((self.vk - vk) / self.vk)**2) / len(self.k)

