import logging
import os
import shutil

import structopt
import structopt.tools.lammps


def run(parameters, individual, relax, use_mpi4py):
    logger = logging.getLogger('by-rank')
    structure = individual #structure = compose_structure(individual)  # TODO
    # Currently turning keep_tmp_files on breaks lammps. Until then,
    # manually delete
    keep_files = parameters.setdefault('keep_files', True)
    temp_parameters = parameters.copy()
    temp_parameters['keep_files'] = True

    calculator = setup_lammps(parameters=temp_parameters, relax=True, use_mpi4py=use_mpi4py)
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
        logger.info('Finished relaxation of individual{0} @ rank {1}: energy = {2}'.format(individual.index, logging.parameters.rank, energy))
    except Exception as error:
        logger.critical('Error in energy evaluation: {0}'.format(error), exc_info=True)
        # Copy files to TroubledLammps directory
        path = os.path.join(cwd,'TroubledLammps')
        os.makedirs(path, exist_ok=True)
        calc = structure.get_calculator()
        for f in [calc.lammps_trj, calc.lammps_in, calc.lammps_log, calc.lammps_data]:
            try:
                shutil.copyfile(f, os.path.join(path, os.path.basename(f)))
            except FileNotFoundError:
                continue
        raise RuntimeError from error

    #individual, buli = decompose_structure(new_structure, individual)  # TODO
    buli = None
    individual = structure

    #individual.buli = buli
    #individual.energy = energy
    individual.LAMMPS = energy
    #individual.pressure = pressure
    #individual.volume = volume

    #structure.calc.clean()
    if not keep_files:
        for f in [calculator.lammps_trj,
                  calculator.lammps_in,
                  calculator.lammps_log,
                  calculator.lammps_data]:
            try:
                os.remove(f)
            except FileNotFoundError:
                continue

    if relax:
        return individual
    else:
        return energy


def setup_lammps(parameters, relax, use_mpi4py):
    logger = logging.getLogger('by-rank')

    # We need to get the ordered list of symbols from the potential file
    potential = structopt.io.eam.read_eam(parameters["potential_file"], kind=parameters["pair_style"])
    ordered_symbols = potential[1][0]  # https://github.com/libAtoms/matscipy/blob/master/matscipy/calculators/eam/io.py

    # TODO I only copied one of the pair_style options
    if parameters["pair_style"] == 'eam/alloy':
        pair_coeff = '* * {0} {1}'.format(parameters["potential_file"], ' '.join(ordered_symbols))
        masses = potential[1][2]
        masses = ["{i} {m}".format(i=i, m=m) for i, m in enumerate(masses, start=1)]
        lammps_parameters = {
            'pair_style': parameters["pair_style"],
            'pair_coeff': [pair_coeff],
            'mass': masses
        }
        files = [parameters["potential_file"]]
    elif parameters["pair_style"] == 'eam':
        pair_coeff = '* * {0}'.format(parameters["potential_file"])
        lammps_parameters = {
            'pair_style': parameters["pair_style"],
            'pair_coeff': [pair_coeff],
        }
        files = [parameters["potential_file"]]

    if parameters["minimize"] != None:
        try:
            lammps_parameters['mass'][len(lammps_parameters['mass'])-1] += '\nmin_style {0}'.format(parameters["min_style"])
        except KeyError:
            lammps_parameters['pair_coeff'][-1] += '\nmin_style {0}'.format(parameters["min_style"])
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
        #path = os.path.join(os.getcwd(), '{0}-rank0'.format(logging.parameters.output_filename))
        path = os.path.join(os.getcwd(), logging.parameters.output_filename)
        rank = logging.parameters.rank
        logger.debug('Setting up directory for keeping LAMMPS files')
        os.makedirs(os.path.join(path, 'LAMMPSFiles'), exist_ok=True)

        # Update kwargs
        kwargs['keep_tmp_files'] = True
        if use_mpi4py:
            tmp_dir = os.path.join(os.path.join(path, 'LAMMPSFiles'), 'rank-{0}'.format(rank))
        else:
            tmp_dir = os.path.join(path, 'LAMMPSFiles')
        kwargs['tmp_dir'] = tmp_dir

    return structopt.tools.lammps.LAMMPS(**kwargs)

