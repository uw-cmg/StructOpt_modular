"""This contains a StructOpt calculator for submitting jobs to queue, tracking their progress, and reading their output. Currently only works with genetic.py"""

import os
import re
import json
from operator import itemgetter

import numpy as np

class StructOpt(object):

    def __init__(self, calcdir=None, optimizer=None, parameters=None, submit_parameters=None):

        if calcdir is None:
            self.calcdir = os.getcwd()
        else:
            self.calcdir = os.path.expandvars(calcdir)
        self.abs_calcdir = os.path.abspath(self.calcdir)
        self.cwd = os.getcwd()
        self.optimizer = optimizer
        self.parameters = parameters
        self.submit_parameters = submit_parameters

        # Set default analysis parameters
        self.fitness = None
        self.log_dirs = None
        self.log_dir = None

        self.system_name = os.path.basename(self.calcdir)

    def __enter__(self):
        """On enter, make sure directory exists. Create it if necessary
        and change into the directory. Then return the calculator."""

        # Make directory if it doesn't already exist
        if not os.path.isdir(self.calcdir):
            os.makedirs(self.calcdir)

        # Now change into the new working directory
        os.chdir(self.calcdir)
        self.initialize()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """On exit, change back to the original directory"""

        os.chdir(self.cwd)

        return

    def initialize(self):
        """Gets the status of the calculation. Meant to be run when within
        the calculation directory"""

        pass

    def run(self):
        """Runs the job as is in the current directory."""

        pass

    def read_runs(self):
        """Stores the output directory in the order in which they were run"""

        log_dirs = [d for d in os.listdir('.') if os.path.isdir('./' + d)]
        pattern = r'logs(.*)'
        log_times = [int(re.match(pattern, d, re.I|re.M).group(1)) for d in log_dirs
                     if re.match(pattern, d, re.I|re.M)]

        log_dirs_times = zip(log_dirs, log_times)
        log_dirs_times = sorted(log_dirs_times, key=itemgetter(1))
        log_dirs, log_times = zip(*log_dirs_times)

        self.log_dirs = log_dirs

        return

    def set_run(self, run_number):
        """Sets the get and read functions on a certain run number.

        Parameters
        ----------
        run_number : int
            The run number we wish to extra data out of. Normal indexing rules
            apply, so run_number = -1 is the last (most recent) run, while 0
            is the first (earliest) run.
        """

        if self.log_dirs is None:
            self.read_runs()

        new_log_dir = self.log_dirs[run_number]
        if new_log_dir is not self.log_dir:
            self.clear_data()

        self.log_dir = new_log_dir

    def clear_data(self):
        """Run to clear the stored data."""

        self.fitnesses = None

        return

    def submit(self):
        """Submits the job to the queue"""

        self.write_input()
        self.write_submit()

        

        return

    def write_input(self):
        """Writes the parameters input to a json file"""

        with open('structopt.in.json', 'w') as f:
            json.dump(self.parameters, f, indent=4, sort_keys=True)

        return

    def write_submit(self):
        """This function writes the submit script."""

        from .rc import QUEUE_OPTIONS as queue_options
        from .rc import RUN_OPTIONS as run_options        
        from .rc import CUSTOM_LINES as custom_lines

        # Gather variables from rc file
        queue_system = 'PBS'
        options = queue_options[queue_system]
        submit = self.submit_parameters

        prefix = options['prefix']
        job_name = options['job_name'].format(submit['job_name'])
        queue = options['queue'].format(submit['queue'])
        nodes_cores = options['nodes_cores'].format(submit['nodes'], submit['cores'])
        total_cores = submit['nodes'] * submit['cores']
        walltime = options['walltime'].format(submit['walltime'])
        misc = options['misc']

        mpirun = run_options['mpirun']
        python = run_options['python']
        optimizer = os.path.expandvars('$STRUCTOPT_HOME/structopt/{}.py'.format(self.optimizer))
        input_file = 'structopt.in.json'

        # Write the submit script
        script = '#!/bin/bash\n\n'
        
        script += '{prefix} {job_name}\n'.format(**locals())
        script += '{prefix} {queue}\n'.format(**locals())
        script += '{prefix} {nodes_cores}\n'.format(**locals())
        script += '{prefix} {walltime}\n'.format(**locals())
        script += '{prefix} {misc}\n\n'.format(**locals())

        script += '{custom_lines}\n\n'.format(**locals())

        script += 'cd {}\n\n'.format(self.abs_calcdir)

        script += '{mpirun} -n {total_cores} {python} {optimizer} {input_file}'.format(**locals())

        with open('submit.sh', 'w') as f:
            f.write(script)

        return

    def read_fitness(self):
        """Reads fitness.log and stores the data"""

        all_fitness = []
        pattern = '.* Generation (.*), Individual (.*): (.*)'

        current_population = []
        current_generation = 0
        with open('{}/fitnesses.log'.format(self.log_dir)) as fitness_file:
            for line in fitness_file:

                # Get the data from the line
                match = re.match(pattern, line, re.I|re.M)
                if match:
                    generation = int(match.group(1))
                    individual = int(match.group(2))
                    fitness = float(match.group(3))
                else:
                    continue

                # If we get here, we've finished reading all fitnesses
                # for a current generation
                if generation > current_generation:
                    current_generation = generation
                    all_fitness.append(current_population)
                    current_population = []

                # Make sure the fitness list is big enough
                if len(current_population) < individual + 1:
                    current_population += [None]*(individual + 1 - len(current_population))

                current_population[individual] = fitness

        self.fitness = np.array(all_fitness)

        return

    def get_fitnesses(self):
        """Returns a list of fitnesses of all individuals in all generations.

        Output
        ------
        out : N x M numpy.ndarray
            N = Number of generations
            M = Number of individuals in each generation
            out[i][j] gives fitness of individual j in generation i.
        """

        if self.fitness is None:
            self.read_fitness()

        return self.fitness

    def get_avg_fitnesses(self):
        """Returns a list of the average fitness of each generation

        Output
        ------
        out : numpy.ndarray
            N = Number of generations
            out[i] gives average fitness of generation i
        """

        if self.fitness is None:
            self.read_fitness()

        return np.average(self.fitness, axis=1)

    def get_min_fitnesses(self):
        """Returns a list of the minimum fitness of each generation

        Output
        ------
        out : N sized numpy.ndarray
            N = Number of generations
            out[i] gives  fitness of generation i
        """

        if self.fitness is None:
            self.read_fitness()

        return np.amin(self.fitness, axis=1)

    def get_max_fitnesses(self):
        """Returns a list of the minimum fitness of each generation

        Output
        ------
        out : N sized numpy.ndarray
            N = Number of generations
            out[i] gives fitness of generation i
        """

        if self.fitness is None:
            self.read_fitness()

        return np.amax(self.fitness, axis=1)

    def get_stdev_fitness(self):
        """Returns the fitness standard deviation of each generation

        Output
        ------
        out : N sized numpy.ndarray
            N = Number of generations
            out[i] gives standard devaiation of generation i
        """

        if self.fitness is None:
            self.read_fitness()

        return np.std(self.fitness, axis=1)
