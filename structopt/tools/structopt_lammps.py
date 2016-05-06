import logging
import os
import shutil

import structopt
import structopt.tools.lammps


def run(parameters, individual, relax):
    logger = logging.getLogger('by-rank')
    structure = individual #structure = compose_structure(individual)  # TODO
    calculator = setup_lammps(parameters, relax=True)
    structure.set_calculator(calculator)
    structure.set_pbc(True)

    # Perform Energy Minimization
    try:
        cwd = os.getcwd()
        # Run LAMMPS
        result = structure.calc.calculate(structure)
        new_structure = result['atoms']  # TODO This is not a full Individual object, it's just an ASE atoms object I think; need to modify the input structure I think?
        new_structure.set_pbc(True)
        pea = result['pea']
        energy = result['thermo'][-1]['pe']
        pressure = 0  # should be modified if enthalpy_fit  TODO I don't know what this means (bc of divide by zero error maybe?)
        volume = new_structure.get_volume()
        logger.info('Finished relaxation of individual{0} @ rank {1}: energy = {2}'.format(individual.index, structopt.parameters.globals.rank, energy))
    except Exception as error:
        logger.critical('Error in energy evaluation: {0}'.format(error), exc_info=True)
        # Copy files to TroubledLammps directory
        path = os.path.join(cwd,'TroubledLammps')
        if not os.path.exists(path):
            os.mkdir(path)
        calc = structure.get_calculator()
        shutil.copyfile(calc.lammps_traj, os.path.join(path, os.path.basename(calc.lammps_traj)))
        shutil.copyfile(calc.lammps_in, os.path.join(path, os.path.basename(calc.lammps_in)))
        shutil.copyfile(calc.lammps_log, os.path.join(path, os.path.basename(calc.lammps_log)))
        shutil.copyfile(calc.lammps_data, os.path.join(path, os.path.basename(calc.lammps_data)))
        raise RuntimeError from error

    #individual, buli = decompose_structure(new_structure, individual)  # TODO
    buli = None
    individual = structure

    #individual.buli = buli
    #individual.energy = energy
    individual.LAMMPS = energy
    #individual.pressure = pressure
    #individual.volume = volume

    structure.calc.clean()

    if relax:
        return individual
    else:
        return energy


def setup_lammps(parameters, relax):
    logger = logging.getLogger('by-rank')
    # TODO I only copied one of the pair_style options
    if parameters["pair_style"] == 'eam/alloy':
        parcoff = '* * {0}'.format(parameters["pot_file"])
        for one in atomlist:
            parcoff += ' {0}'.format(one[0])
        pair_coeff = [parcoff]
        lammps_parameters = {'pair_style': parameters["pair_style"],
                      'pair_coeff': pair_coeff}
        files = [parameters["potential_file"]]
    elif parameters["pair_style"] == 'eam':
        pair_coeff = ['* * {0}'.format(parameters["potential_file"])]
        lammps_parameters = {'pair_style': parameters["pair_style"],
                      'pair_coeff': pair_coeff}
        files = [parameters["potential_file"]]

    if parameters["minimize"] != None:
        try:
            lammps_parameters['mass'][len(lammps_parameters['mass'])-1] += '\nmin_style {0}'.format(parameters["min_style"])
        except KeyError:
            lammps_parameters['pair_coeff'][0] += '\nmin_style {0}'.format(parameters["min_style"])
        lammps_parameters['minimize'] = parameters["minimize"]

    if not relax:
        lammps_parameters['minimize'] = "1e-8 1e-8 0 0"  # Set maximim steps to 0
    lammps_parameters['thermosteps'] = parameters["thermo_steps"]
    
    # Set up kwargs for LAMMPS calculator
    kwargs = {
            'parameters': lammps_parameters
    }

    if files != None:
        kwargs['files'] = files

    if parameters["keep_files"]:
        # Set up directory for saving files
        #path = os.path.join(os.getcwd(), '{0}-rank0'.format(structopt.parameters.globals.output_filename))
        path = os.path.join(os.getcwd(), structopt.parameters.globals.output_filename)
        rank = structopt.parameters.globals.rank
        logger.debug('Setting up directory for keeping LAMMPS files')
        if structopt.parameters.globals.USE_MPI4PY:
            if not os.path.exists(os.path.join(path, 'LAMMPSFiles')):
                os.mkdir(os.path.join(path, 'LAMMPSFiles'))
                logger.info('Making directory: {0}'.format(os.path.join(path, 'LAMMPSFiles')))

        # Update kwargs
        kwargs['keep_tmp_files'] = True
        if structopt.parameters.globals.USE_MPI4PY:
            tmp_dir = os.path.join(os.path.join(path, 'LAMMPSFiles'), 'rank-{0}'.format(rank))
        else:
            tmp_dir = os.path.join(path, 'LAMMPSFiles')
        kwargs['tmp_dir'] = tmp_dir

    return structopt.tools.lammps.LAMMPS(**kwargs)

