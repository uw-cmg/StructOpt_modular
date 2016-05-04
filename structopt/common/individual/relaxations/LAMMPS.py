import structopt
from StructOpt.tools.lammps import LAMMPS


class LAMMPS(object):
    def __init__(self):
        self.parameters = structopt.parameters.relaxations.LAMMPS


    def get_command(self, individual):
        raise NotImplementedError


    def relax(self, individual):
        logger = logging.getLogger('by-rank')
        structure = compose_structure(individual)  # TODO
        calculator = self.setup_lammps(self.parameters, relax)  # TODO
        solution.set_calculator(calculator)
        solution.set_pbc(True)

        # Perform Energy Minimization
        try:
            cwd = os.getcwd()
            # Run LAMMPS
            result = structure.calc.calculate(structure.copy())
            new_structure = result['atoms']  # TODO This is not a full Individual object, it's just an ASE atoms object I think; need to modify the input structure I think?
            new_structure.set_pbc(True)
            pea = result['pea']
            energy = result['thermo'][-1]['pe']
            pressure = 0  # should be modified if enthalpy_fit  TODO I don't know what this means (bc of divide by zero error maybe?)
            volume = new_structure.get_volume()
            logger.info('Finished relaxation of individual{0} @ rank {1}: energy = {2}'.format(individual.index, rank, energy))
        except Exception as error:
            logger.critical('Error in energy evaluation: {0}'.format(error), exc_info=True)
            # Copy files to TroubledLammps directory
            path = os.path.join(cwd,'TroubledLammps')
            if not os.path.exists(path):
                os.mkdir(path)
            shutil.copyfile(calc.trajfile, os.path.join(path, os.path.basename(calc.trajfile)))
            shutil.copyfile(calc.infile, os.path.join(path, os.path.basename(calc.infile)))
            shutil.copyfile(calc.logfile, os.path.join(path, os.path.basename(calc.logfile)))
            shutil.copyfile(calc.datafile, os.path.join(path, os.path.basename(calc.datafile)))
            raise RuntimeError('{0}:{1}'.format(Exception, error)) from error

        individual, buli = decompose_structure(new_structure, individual)  # TODO
        individual.buli = buli
        individual.energy = energy
        individual.pressure = pressure
        individual.volume = volume

        calc.clean()

        if relax:
            return individual
        else:
            return energy


    def setup_lammps(parameters, relax):
        # TODO there are a bunch of bugs/stuff I don't understand in here
        # TOD I only copied one of the pair_style options
        if parameters["pair_style"] == 'eam/alloy':
            parcoff = '* * {0}'.format(parameters["pot_file"])
            for one in atomlist:
                parcoff += ' {0}'.format(one[0])
            pair_coeff = [parcoff]
            lammps_parameters = {'pair_style': parameters["pair_style"],
                          'pair_coeff': pair_coeff}
            files = [parameters["pot_file"]]

        if parameters["minimize"] != None:
            try:
                lammps_parameters['mass'][len(lammps_parameters['mass'])-1] += '\nmin_style {0}'.format(parameters["min_style"])
            except KeyError:
                lammps_parameters['pair_coeff'][0] += '\nmin_style {0}'.format(parameters["min_style"])
            lammps_parameters['minimize'] = parameters["minimize"]

        if not relax:
            lammps_parameters['minimize'] = "1e-8 1e-8 0 0"
        lammps_parameters['thermosteps'] = parameters["thermo_steps"]
        
        # if parameters["keep_files"]:
        # TODO I did not copy all this

        return LAMMPS(parameters=lammps_parameters, files=files)

