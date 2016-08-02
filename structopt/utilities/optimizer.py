"""This contains a StructOpt calculator for submitting jobs to queue, tracking their progress, and reading their output."""

import os
import re
import json
from copy import deepcopy
from operator import itemgetter
import subprocess
from subprocess import PIPE

import numpy as np
import ase

from ..common.individual import Individual
from .exceptions import StructOptUnknownState, StructOptRunning, StructOptQueued, StructOptSubmitted

class StructOpt(object):

    def __init__(self, calcdir=None, optimizer=None, parameters=None,
                 submit_parameters=None):

        # Initialize inputs
        if calcdir is None:
            self.calcdir = os.getcwd()
        else:
            self.calcdir = os.path.expandvars(calcdir)
        self.optimizer = optimizer
        self.parameters = deepcopy(parameters)
        self.submit_parameters = deepcopy(submit_parameters)
        if 'job_name' not in self.submit_parameters:
            self.submit_parameters['job_name'] = self.calcdir

        self.path = os.path.abspath(self.calcdir)
        self.cwd = os.getcwd()
        self.system_name = os.path.basename(self.calcdir)

        # Initialize results dictionaries
        self.fitness = None
        self.genenerations = None
        self.individuals = None
        self.log_dirs = None
        self.log_dir = None

        # If totally clean calculation
        if not os.path.isdir(self.path):
            self.status = 'clean'
            os.makedirs(self.path)

        # If directory exists
        elif not os.path.isfile(os.path.join(self.path, 'structopt.in.json')):
            self.status = 'clean'

        # If input exists but was never run
        elif (os.path.isfile(os.path.join(self.path, 'structopt.in.json'))
              and not self.job_in_queue(os.path.join(self.path, 'jobid'))
              and not self.read_runs()):
            self.read_input()
            self.status = 'initialized'

        # If the job is running
        elif (os.path.isfile(os.path.join(self.path, 'structopt.in.json'))
              and self.job_in_queue(os.path.join(self.path, 'jobid'))):
            self.read_input()
            self.status = 'running'

        # If the job is done, check the output
        elif (os.path.isfile(os.path.join(self.path, 'structopt.in.json'))
              and not self.job_in_queue(os.path.join(self.path, 'jobid'))
              and self.read_runs()):
            self.read_input()
            self.status = self.check_run() # Can be "done" or "error"
            self.read_generations()

    def restart(self): # TODO
        """Loads up the last generation of a previous run and modifies 
        the self.parameters to load up those structures on the next run"""

        pass

    def job_in_queue(self, jobid='jobid'):
        '''return True or False if the directory has a job in the queue'''
        if not os.path.exists(jobid):
            return False

        with open(os.path.join(self.path, 'jobid')) as f:
            jobid = f.readline().split()[-1]

        # Behavior will depend on whether we are in the slurm or pbs environment
        if self.submit_parameters['system'] == 'PBS':
            jobids_in_queue = subprocess.check_output('qselect')
            jobids_in_queue = [job.decode('utf-8') for job in jobids_in_queue.split()]
        else:
            raise NotImplemented(self.submit_parameters['system'], 'not implemented yet')

        if jobid in jobids_in_queue:
            # get details on specific jobid
            output, error = subprocess.Popen(['qstat', format(jobid)], stdout=PIPE).communicate()
            if error is None:
                fields = output.decode('utf-8').split('\n')[-2].split()
                job_status = fields[4]
                if job_status == 'C':
                    return False
                else:
                    return True
            return False
        else:
            return False


    ####################################################################
    ### Calculation methods. Includes write, run, and submit scripts ###
    ####################################################################

    def optimize(self, run_method='submit', rerun=False, restart=False):
        """Runs the optimizer"""

        run_method = getattr(self, run_method)

        if restart:
            self.restart()

        if self.status == 'clean' or rerun:
            run_method()

        return

    def run(self): # TODO
        """Runs the job as is in the current directory."""

        self.write_input()
        self.write_submit()

        return

    def submit(self):
        """Submits the job to the queue. Do this in the calculation
        directory"""

        from .rc import QUEUE_OPTIONS as queue_options

        self.write_input()
        self.write_submit()

        submit_cmd = queue_options[self.submit_parameters['system']]['submit']

        os.chdir(self.path)

        p = subprocess.Popen([submit_cmd, 'submit.sh'], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()

        with open('jobid', 'wb') as f:
            f.write(out)

        os.chdir(self.cwd)

        raise StructOptSubmitted(out)

    def write_input(self):
        """Writes the parameters input to a json file"""

        input_file = os.path.join(self.path, 'structopt.in.json')
        with open(input_file, 'w') as f:
            json.dump(self.parameters, f, indent=4, sort_keys=True)

        return

    def write_submit(self):
        """This function writes the submit script."""

        from .rc import QUEUE_OPTIONS as queue_options
        from .rc import RUN_OPTIONS as run_options
        from .rc import CUSTOM_LINES as custom_lines

        # Gather variables from rc file
        submit = self.submit_parameters
        queue_system = submit['system']
        options = queue_options[queue_system]

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

        script += 'cd {}\n\n'.format(self.path)

        if total_cores > 1:
            script += '{mpirun} -n {total_cores} {python} {optimizer} {input_file}'.format(**locals())
        else:
            script += '{python} {optimizer} {input_file}'.format(**locals())

        submit_file = os.path.join(self.path, 'submit.sh')
        with open(submit_file, 'w') as f:
            f.write(script)

        return

    #############################################################
    ### Postprocessing methods. Includes read and get scripts ###
    #############################################################

    def read_runs(self):
        """Stores the output directory in the order in which they were run. 
        Sets StructOpt to read results form most recent by default"""

        dirs = [d for d in os.listdir(self.path)
                if os.path.isdir(os.path.join(self.path, d))]
        pattern = r'logs(.*)'
        log_times = [int(re.match(pattern, d, re.I|re.M).group(1)) for d in dirs
                     if re.match(pattern, d, re.I|re.M)]
        log_dirs = [os.path.join(self.path, 'logs{}'.format(t)) for t in log_times]

        if len(log_dirs) > 0:
            log_dirs_times = zip(log_dirs, log_times)
            log_dirs_times = sorted(log_dirs_times, key=itemgetter(1))
            log_dirs, log_times = zip(*log_dirs_times)
            self.log_dirs = log_dirs
            self.set_run(-1)
            return True
        else:
            self.log_dirs = []
            return False

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

    def check_run(self):
        """Check the stdout to see if a run is complete. NOTE, only works
        for jobs submitted to the queue and jobs where the stdout is
        saved to stdout.txt"""

        if 'stdout.txt' in os.listdir(self.path):
            out_file = os.path.join(self.path, 'stdout.txt')
        else:
            files = os.listdir(self.path)
            pattern = '.*.o(.*)'
            files = [f for f in files if re.match(pattern, f, re.I|re.M)]
            files = [f for f in files if re.match(pattern, f, re.I|re.M).group(1).isnumeric()]
            jobids = [int(re.match(pattern, f, re.I|re.M).group(1)) for f in files]
            files_jobids = zip(files, jobids)
            files_jobids = sorted(files_jobids, key=lambda i: i[1])
            out_file = files_jobids[-1][0]

        with open(os.path.join(self.path, out_file)) as f:
            last = None
            for last in (line for line in f if line.rstrip('\n')):
                pass

        # Check if the job ran error free. Happens with "Finished!" gets
        # printed and if the job was killed due to walltime
        if ('Finished!' in last):
            return 'done'
        elif ('walltime' in last):
            return 'timeout'
        else:
            return 'error'

    def read_generations(self):
        """Determines which generations exists and initializes
        a list to store the individuals"""

        # Get list of generations available to be read for output
        dirs = [d for d in os.listdir(os.path.join(self.log_dir, 'XYZs'))]
        pattern = r'generation(.*)'
        generations = [int(re.match(pattern, d, re.I|re.M).group(1)) for d in dirs]
        self.generations = sorted(generations)

        # Initialize dictionary for storing XYZ coordinates
        self.populations = [None for generation in range(max(generations) + 1)]

        return

    def read_population(self, generation):
        """Reads and stores extended .xyz files into self.individuals dictionary"""

        # Get a list of individuals available in the generation
        generation_dir = os.path.join(self.log_dir,
                                      'XYZs/generation{}'.format(generation))

        population = []
        for f in os.listdir(generation_dir):
            pattern = r'individual(.*).xyz'
            index = int(re.match(pattern, f, re.I|re.M).group(1))
            atoms = ase.io.read(os.path.join(generation_dir, f))

            individual = Individual(index=index)
            individual.extend(atoms)
            individual.set_pbc(atoms.get_pbc())
            individual.set_cell(atoms.get_cell())
            population.append(individual)

            # TODO: Probably need to add methods for reading other properties

        population.sort(key=lambda individual: individual.index)
        self.populations[generation] = population

        return

    def get_population(self, generation=-1):
        """Returns a list of Individual objects, ordered by their index,
        for a given generation. Note this is NOT a structopt Population.

        Parameters
        ----------
        generation : int
            Specifies which generation to return.
        """

        if generation < 0:
            generation = max(self.generations) + generation

        if generation not in self.generations:
            raise IOError('Generation {} not found'.format(generation))

        self.read_population(generation)
        return self.populations[generation]

    def get_all_populations(self):
        """Returns all of populations stored"""

        for generation in self.generations:
            self.read_population(generation)

        return self.populations

    def clear_data(self): # TODO
        """Run to clear the stored data."""

        self.fitnesses = None

        return

    def check_for_errors(self): # TODO
        """After obtaining list of log directories, eliminate directories
        that returned an error or is not done"""

        pass

    def read_input(self): # TODO
        """Read the input and store them self.parameter"""

        pass

    def read_fitness(self):
        """Reads fitness.log and stores the data"""

        all_fitnesses = {'total': []}
        current_fitnesses = {'total': []}
        modules = self.parameters['fitnesses']

        pattern = '.* Generation (.*), Individual (.*):'
        for module in modules:
            pattern += ' (.*): (.*)'
            all_fitnesses.update({module: []})
            current_fitnesses.update({module: []})

        current_generation = 0
        with open(os.path.join(self.log_dir, 'fitnesses.log')) as fitness_file:
            for line in fitness_file:

                # Get the data from the line
                match = re.match(pattern, line, re.I|re.M)
                if match:
                    generation = int(match.group(1))
                    index = int(match.group(2))
                    fitness = {match.group(2*i + 3): match.group(2*i + 4) for i in range(len(modules))}
                    for module in fitness:
                        try:
                            fitness[module] = float(fitness[module])
                        except ValueError:
                            fitness[module] = np.nan
                else:
                    continue

                # If we get here, we've finished reading all fitnesses
                # for a current generation
                if generation > current_generation:
                    current_generation = generation
                    all_fitnesses['total'].append(current_fitnesses['total'])
                    for module in modules:
                        all_fitnesses[module].append(current_fitnesses[module])
                    current_fitnesses = {module: [] for module in current_fitnesses}

                total_fit = 0
                # Store data of current individual
                for module in modules:
                    # Make sure the fitness list is big enough
                    if len(current_fitnesses[module]) < index + 1:
                        current_fitnesses[module].extend([None]*(index + 1 - len(current_fitnesses[module])))
                    current_fitnesses[module][index] = fitness[module]
                    total_fit += fitness[module]

                if len(current_fitnesses['total']) < index + 1:
                    current_fitnesses['total'].extend([None]*(index + 1 - len(current_fitnesses['total'])))
                current_fitnesses['total'][index] = total_fit

            # Append the last generation
            all_fitnesses['total'].append(current_fitnesses['total'])
            for module in modules:
                all_fitnesses[module].append(current_fitnesses[module])

        self.fitness = {module: np.array(all_fitnesses[module]) for module in all_fitnesses}

        return

    def get_fitnesses(self, module='all'):
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

        if module == 'all':
            return self.fitness

        return self.fitness[module]

    def get_avg_fitnesses(self, module='all'):
        """Returns a list of the average fitness of each generation

        Output
        ------
        out : numpy.ndarray
            N = Number of generations
            out[i] gives average fitness of generation i
        """

        if self.fitness is None:
            self.read_fitness()

        if module == 'all':
            self.fitness = {module: np.average(self.fitness[module], axis=1) for module in self.fitness}

        return np.average(self.fitness[module], axis=1)

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

        if module == 'all':
            self.fitness = {module: np.amin(self.fitness[module], axis=1) for module in self.fitness}

        return np.amin(self.fitness[module], axis=1)            

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

        if module == 'all':
            self.fitness = {module: np.amax(self.fitness[module], axis=1) for module in self.fitness}

        return np.amax(self.fitness[module], axis=1)            


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

        if module == 'all':
            self.fitness = {module: np.std(self.fitness[module], axis=1) for module in self.fitness}

        return np.std(self.fitness[module], axis=1)
