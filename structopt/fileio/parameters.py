import json
import logging
import os
from structopt.tools.dictionaryobject import DictionaryObject


def read(input):
    """Sets StructOpt parameters from a dictionary or filename"""

    # Read and store the input
    if isinstance(input, dict):
        parameters = input
    elif isinstance(input, str) and os.path.exists(input):
        parameters = json.load(open(input))
    else:
        raise IOError('Error in input or input file')

    parameters.globals.setdefault('USE_MPI4PY', True)
    parameters.globals.setdefault('output_filename', 'Output')

    return parameters


def write(parameters):
    output = logging.getLogger('output')
    output.info('Current parameters:')
    ouptut.info(parameters.to_json(sorted=True, indent=2))

    return

def set_default(parameters):
    logger = logging.getLogger('default')

    parameters.globals.setdefault('parallel', False)

    if 'fitnesses' not in parameters.globals or not parameters.globals['fitnesses']:
        raise ValueError('Fitnesses must be specified in the parameter file.')

	parameters.globals.setdefault('weights', [1.0 for _ in parameters.globals['fitnesses']])

    if 'relaxations' not in parameters.globals or not parameters.globals['relaxations']:
        raise ValueError('Relaxations must be specified in the parameter file.')

    if 'structure' not in parameters:
        logger.critical("Input file/dictionary must include a structure for the simulation as 'structure':'Cluster/Crystal/Defect'")
        logger.debug("Current parameters include:\n"+repr(parameters))
        raise RuntimeError("Input file/dictionary must include a structure for the simulation as 'structure':'Cluster/Crystal/Defect'")

    if 'atomlist' not in parameters:
        logger.critical("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx', Concentration, Mass, Chemical Potential)]")
        logger.debug("Current parameters include:\n" + repr(parameters))
        raise RuntimeError("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx', Concentration, Mass, Chemical Potential)]")
    else:
        if isinstance(parameters['atomlist'], str):
            parameters['atomlist'] = eval(parameters['atomlist'].strip())
        if not isinstance(parameters['atomlist'], list):
            logger.critical('Something is wrong with atomlist parameter: {0}'.format(parameters['atomlist']))
            raise RuntimeError("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx', Concentration, Mass, Chemical Potential)]: {0}".format(parameters['atomlist']))

    for one in parameters['atomlist']:
        if len(one) != 4:
            # Stop program if atom list parameter not properly formatted
            logger.critical('Format of atom list not correct. Must be [(symbol, concentration, mass, potential)]')
            logger.debug('Issue in section : {0}'.format(one))
            logger.debug('Current atomlist is formatted as : {0}'.format(parameters['atomlist']))
            raise RuntimeError('Format of atom list not correct. Must be [(symbol, concentration, mass, potential)]')

    if 'natoms' not in parameters:
        parameters['natoms'] = int(sum([abs(c) for ind, c, m, u in parameters['atomlist']]))
        logger.warning('Number of atoms in simulation not set')
        logger.warning('Assuming natoms = {0}'.format(parameters['natoms']))

    parameters['atomlist'] = check_atomlist_concentration(parameters['atomlist'], parameters['natoms'])

    if 'optimizer_type' not in parameters:
        logger.info('optimizer_type not set.  Default values set to GA.')
    parameters.setdefault('optimizer_type', 'GA')

    if parameters['optimizer_type'] == 'GA':
        nindiv = 10
        genealogy = True
        nbests = 100
        algtype = 'lambda+mu'
        cxpb = 0.8
        mutpb = 0.15
        natselectscheme = 'tournament'
    else:
        nindiv = 1
        genealogy = False
        nbests = 100
        algtype = parameters['optimizer_type']
        cxpb = 0.0
        mutpb = 1.0
    if parameters['optimizer_type'] == 'SA':
        natselectscheme = 'metropolis'
        predator = 'adapting'
    elif parameters['optimizer_type'] == 'BH':
        natselectscheme = 'metropolis'
    else:
        natselectscheme = 'best'
    if 'nindiv' not in parameters:
        parameters['nindiv'] = nindiv
        logger.info('Setting number of individuals in population (nindiv) = {0}'.format(parameters['nindiv']))

    # Parameters for output
    if 'genealogy' not in parameters:
        parameters['genealogy'] = genealogy
        logger.info('Setting genealogy = {0}'.format(parameters['genealogy']))
    if 'output_format' not in parameters:
        parameters['output_format'] = 'fitness'
        logger.info('Setting output format = {0}'.format(parameters['output_format']))
    if 'allenergyfile' not in parameters:
        parameters['allenergyfile'] = False
        logger.info('Setting allenergyfile = {0}'.format(parameters['allenergyfile']))
    if 'best_inds_list' not in parameters:
        parameters['best_inds_list'] = True
        logger.info('Setting best_inds_list = {0}'.format(parameters['best_inds_list']))
    if 'number_of_bests' not in parameters:
        parameters['number_of_bests'] = nbests
        logger.info('Setting number_of_bests = {0}'.format(parameters['number_of_bests']))
    if 'indiv_defect_write' not in parameters:
        parameters['indiv_defect_write'] = False
        logger.info('Setting indiv_defect_write = {0}'.format(parameters['indiv_defect_write']))
    if 'vacancy_output' not in parameters:
        parameters['vacancy_output'] = False
        logger.info('Setting vacancy_output = {0}'.format(parameters['vacancy_output']))
    if 'restart_optimizer' not in parameters:
        parameters['restart_optimizer'] = False
    if 'restart_files' not in parameters:
        parameters['restart_files'] = True
        logger.info('Setting restart_files = {0}'.format(parameters['restart_files']))
    if 'indiv_write' not in parameters:
        parameters['indiv_write'] = 'all'
        logger.info('Setting indiv_write = {0}'.format(parameters['indiv_write']))

    # Parameters for post-processing
    if 'lattice_concentration' not in parameters:
        parameters['lattice_concentration'] = False
        logger.info('Setting lattice_concentration = {0}'.format(parameters['lattice_concentration']))
    if 'postprocessing' not in parameters:
        parameters['postprocessing'] = False
        logger.info('Setting postprocessing = {0}'.format(parameters['postprocessing']))
    if 'genealogytree' not in parameters:
        parameters['genealogytree'] = False
        logger.info('Setting genealogytree = {0}'.format(parameters['genealogytree']))

    # Parameters for general algorithm
    if 'seed' not in parameters:
        parameters['seed'] = random.randint(0, 10)
        logger.info('Setting Random number seed (seed) to {0}'.format(parameters['seed']))
    if 'forcing' not in parameters:
        parameters['forcing'] = 'Concentration'
        logger.info('Setting forcing = {0}'.format(parameters['forcing']))
        logger.info('Assuming forcing concentration control')
    if 'debug' not in parameters:
        parameters['debug'] = ['None']
        logger.info('Setting debug = {0}'.format(parameters['debug']))
        if 'None' not in parameters['debug']:
            print '***** DEBUGGING RUN *****'
    if 'algorithm_type' not in parameters:
        parameters['algorithm_type'] = algtype
        logger.info('Setting algorithm type = {0}'.format(parameters['algorithm_type']))
    if 'migration_intervals' not in parameters:
        parameters['migration_intervals'] = 5
        logger.info('Setting migration_intervals = '.format(parameters['migration_intervals']))
    if 'migration_percent' not in parameters:
        parameters['migration_percent'] = 0.05
        logger.info('Setting migration_percent = {0}'.format(parameters['migration_percent']))
    if 'fingerprinting' not in parameters:
        parameters['fingerprinting'] = False
        logger.info('Setting fingerprinting = {0}'.format(parameters['fingerprinting']))
    if 'fpbin' not in parameters:
        parameters['fpbin'] = 0.25
        logger.info('Setting fingerprint bin to {0}'.format(parameters['fpbin']))
    if 'fpcutoff' not in parameters:
        parameters['fpcutoff'] = 15.0
        logger.info('Setting fingerprint cutoff distance to {0}'.format(parameters['fpcutoff']))
    parameters['bulkfp'] = None
    if 'fixed_region' not in parameters:
        parameters['fixed_region'] = False
        logger.info('Setting fixed_region = {0}'.format(parameters['fixed_region']))
    if 'rattle_atoms' not in parameters:
        parameters['rattle_atoms'] = False
        logger.info('Setting rattle_atoms = {0}'.format(parameters['rattle_atoms']))
    if 'constrain_position' not in parameters:
        parameters['constrain_position'] = False
        logger.info('Setting constrain_position = {0}'.format(parameters['constrain_position']))
    if 'restart' not in parameters:
        parameters['restart'] = False
        logger.info('Setting restart = {0}'.format(parameters['restart']))
    if 'restart_ints' not in parameters:
        parameters['restart_ints'] = 0
        if parameters['restart']:
            logger.info('Setting restart_ints = {0}'.format(parameters['restart_ints']))

    # Parameters to generate the population and individual
    if 'r_ab' not in parameters:
        parameters['r_ab'] = 2.5
        logger.info('Setting r_ab = {0}'.format(parameters['r_ab']))
    if 'size' not in parameters:
        parameters['size'] = parameters['natoms']**0.33333*parameters['r_ab']
        logger.info('Setting size to (natoms)^(1/3)*r_ab = {0}'.format(parameters['size']))
    if 'generate_flag' not in parameters:
        parameters['generate_flag'] = 'box'
        logger.info('Setting default generation scheme = {0}'.format(parameters['generate_flag']))
    parameters['solidbulk'] = None
    if 'sf' not in parameters:
        parameters['sf'] = 1.75
        if parameters['structure'] == 'Defect':
            logger.info('Setting size factor for Defect (sf) = {0}'.format(parameters['sf']))
    if 'supercell' not in parameters:
        parameters['supercell'] = (1, 1, 1)
        if parameters['structure'] == 'Defect':
            logger.info('Setting supercell for Defect (supercell) = {0}'.format(parameters['supercell']))
    if 'solidfile' not in parameters:
        if parameters['structure'] == 'Defect':
            logger.critical('Must provide a file for bulk solid if running a defect simulation.')
            raise RuntimeError('Error: Bulk for Defect not specified. Enter name of file for bulk structure as solidfile parameter')
        else:
            parameters['solidfile'] = None
    if 'solidcell' not in parameters:
        if parameters['structure'] == 'Defect':
            logger.warning('Warning cell size for Bulk Solid not specified assuming distance between 1st and 3rd atom')
        parameters['solidcell'] = None
    if 'evalsolid' not in parameters:
        parameters['evalsolid'] = False
        if parameters['structure'] == 'Defect':
            logger.info('Not evaluating Solid')
    if 'finddefects' not in parameters:
        parameters['finddefects'] = True
        if parameters['structure'] == 'Defect':
            logger.info('Setting finddefects = {0}'.format(parameters['finddefects']))
    if 'trackvacs' not in parameters:
        parameters['trackvacs'] = False
        if parameters['structure'] == 'Defect':
            logger.info('Setting trackvacs = {0}'.format(parameters['trackvacs']))
    if 'trackswaps' not in parameters:
        parameters['trackswaps'] = False
        if parameters['structure'] == 'Defect':
            logger.info('Setting trackswaps = {0}'.format(parameters['trackswaps']))
    if 'random_loc_start' not in parameters:
        parameters['random_loc_start'] = False
        if parameters['structure'] == 'Defect':
            logger.info('Setting random_loc_start = {0}'.format(parameters['random_loc_start']))
    if 'random_vac_start' not in parameters:
        parameters['random_vac_start'] = False
        if parameters['structure'] == 'Defect':
            logger.info('Setting random_vac_start = {0}'.format(parameters['random_vac_start']))
    if 'purebulkenpa' not in parameters:
        if parameters['structure'] == 'Defect':
            parameters['purebulkenpa'] = None
    if 'natomsbulk' not in parameters:
        if parameters['structure'] == 'Defect':
            parameters['natomsbulk'] = None
    if 'surfacefile' not in parameters:
        if parameters['structure'] == 'Surface':
            logger.critical('Must provide a file for bulk solid if running a defect simulation.')
            raise RuntimeError('Error: Bulk for Defect not specified. Enter name of file for bulk structure as solidfile parameter')
        else:
            parameters['surfacefile'] = None
    if 'surfacecell' not in parameters:
        parameters['surfacecell'] = None
    if 'surftopthick' not in parameters:
        parameters['surftopthick'] = 0
    if 'cell_shape_options' not in parameters:
        parameters['cell_shape_options'] = ['cubic', 'hexagonal', 'triclinic', \
        'monoclinic', 'orthorhombic', 'tetragonal']
        if parameters['structure'] == 'Crystal':
            logger.info('Assuming following cell shape options: {0}'.format(parameters['cell_shape_options']))
    if 'alloy' not in parameters:
        parameters['alloy'] = True
        logger.info('Setting alloy = {0}'.format(parameters['alloy']))

    if 'large_box_size' not in parameters:
        parameters['large_box_size'] = 500.0
        if parameters['structure'] == 'Cluster':
            logger.info('Setting large_box_size to {0}'.format(parameters['large_box_size']))

    # Parameters for Crossovers
    if 'cxpb' not in parameters:
        parameters['cxpb'] = cxpb
        logger.info('Setting crossover probability (cxpb) = {0}'.format(parameters['cxpb']))
    if 'cx_scheme' not in parameters:
        parameters['cx_scheme'] = 'cxtp'
        logger.info('Assuming two-point crossover.  Setting cx_scheme = {0}'.format(parameters['cx_scheme']))
    if 'selection_scheme' not in parameters:
        parameters['selection_scheme'] = 'tournament2'
        logger.info('Setting selection_scheme = {0}'.format(parameters['selection_scheme']))

    # Parameters for Mutations
    if 'mutpb' not in parameters:
        parameters['mutpb'] = mutpb
        logger.info('Setting mutation probability (mutpb) = {0}'.format(parameters['mutpb']))
    if 'mutation_options' not in parameters:
        if parameters['structure'] == 'Cluster':
            parameters['mutation_options'] = ['lattice_alteration', 'rotation', \
            'permutation', 'scale_size']
        elif parameters['structure'] == 'Crystal':
            parameters['mutation_options'] = ['lattice_alteration', 'rotation', \
            'permutation', 'scale_size', 'cell_shape', 'lammps_box_relax']
        elif parameters['structure'] == 'Defect':
            parameters['mutation_options'] = ['lattice_alteration', 'rotation', 'permutation']
        logger.info('Setting mutations options = {0}'.format(parameters['mutation_options']))
    BHFlag = False
    for one in parameters['mutation_options']:
        if 'basin_hop' in one:
            BHFlag = True
    if 'bh_steps' not in parameters:
        parameters['bh_steps'] = 100
        if BHFlag == True:
            logger.warning('Max steps not specified for Basin Hop mutation, setting bh_steps = {0}'.format(parameters['bh_steps']))
    if 'bh_temp' not in parameters:
        parameters['bh_temp'] = 1000*8.617385692256675e-05
        if BHFlag == True:
            logger.warning('Temperature not set for Basin Hop mutation, setting bh_temp in kT = {0}'.format(parameters['bh_temp']))
    if 'mutant_add' not in parameters:
        parameters['mutant_add'] = False
        logger.info('Setting mutant_add = {0}'.format(parameters['mutant_add']))
    if 'quench_max_temp' not in parameters:
        parameters['quench_max_temp'] = 1000
        if 'quench' in parameters['mutation_options']:
            logger.info('Setting quench_max_temp = {0}'.format(parameters['quench_max_temp']))
    if 'quench_min_temp' not in parameters:
        parameters['quench_min_temp'] = 2
        if 'quench' in parameters['mutation_options']:
            logger.info('Seting quench_min_temp = {0}'.format(parameters['quench_min_temp']))
    if 'quench_step_size' not in parameters:
        parameters['quench_step_size'] = 0.01
        if 'quench' in parameters['mutation_options']:
            logger.info('Setting quench_step_size = {0}'.format(parameters['quench_step_size']))
    if 'quench_n_steps_1' not in parameters:
        parameters['quench_n_steps_1'] = 10000
        if 'quench' in parameters['mutation_options']:
            logger.info('Setting quench_n_steps_1 = {0}'.format(parameters['quench_n_steps_1']))
    if 'quench_n_steps_2' not in parameters:
        parameters['quench_n_steps_2'] = parameters['quench_n_steps_1']*2
        if 'quench' in parameters['mutation_options']:
            logger.info('Setting quench_n_steps_2 = {0}'.format(parameters['quench_n_steps_2']))
    if 'isolate_mutation' not in parameters:
        parameters['isolate_mutation'] = False
        logger.info('Setting isolate_mutation flag = {0}'.format(parameters['isolate_mutation']))

    # Parameters for Selection
    if 'energy_cutoff_factor' not in parameters:
        parameters['energy_cutoff_factor'] = 10.0
        logger.info('Setting energy_cutoff_factor = {0}'.format(parameters['energy_cutoff_factor']))

    if 'constrain_swaps' not in parameters:
        if 'IntSwap' in parameters['mutation_options']:
            parameters['swaplist'] = None
            logger.info('Setting swaplist for IntSwap = None')
        else:
            parameters['swaplist'] = False
            parameters['constrain_swaps'] = False
            logger.info('Setting swaplist = False')
    if 'natural_selection_scheme' not in parameters:
        parameters['natural_selection_scheme'] = natselectscheme
        logger.info('Setting natural_selection_scheme = {0}'.format(parameters['natural_selection_scheme']))
    if 'tournsize' not in parameters:
        parameters['tournsize'] = 3
        if 'tournament' in parameters['selection_scheme'] or 'tournament' in parameters['natural_selection_scheme']:
            logger.info('Setting Tournament size (tournsize) = {0}'.format(parameters['tournsize']))
    if 'fusslimit' not in parameters:
        parameters['fusslimit'] = 10.0
        logger.info('Setting FUSS limit (fusslimit) = {0}'.format(parameters['fusslimit']))
    if 'metropolis_temp' not in parameters:
        parameters['metropolis_temp'] = 30.0
        logger.info('Setting metropolis_temp = {0}'.format(parameters['metropolis_temp']))
    if 'mark' not in parameters:
        parameters['mark'] = None

    # Parameters for Convergence
    if 'convergence_scheme' not in parameters:
        parameters['convergence_scheme'] = 'max_gen'
        logger.info('Setting convergence scheme (convergence_scheme) = {0}'.format(parameters['convergence_scheme']))
    if 'maxgen' not in parameters:
        parameters['maxgen'] = 5
        logger.info('Setting Max Number of generations (maxgen) = {0}'.format(parameters['maxgen']))
    if 'reqrep' not in parameters:
        parameters['reqrep'] = 10
        if 'rep' in parameters['convergence_scheme']:
            logger.info('Setting max number of energy repetitions (reqrep) = {0}'.format(parameters['reqrep']))
    if 'tolerance' not in parameters:
        parameters['tolerance'] = 0.001
        if 'rep' in parameters['convergence_scheme']:
            logger.info('Setting energy tolerance (tolerance) = {0}'.format(parameters['tolerance']))
    if 'predator' not in parameters:
        parameters['predator'] = 'mutation_dups'
        logger.info('Setting predator = {0}'.format(parameters['predator']))
    if 'adaptbegin' not in parameters:
        parameters['adaptbegin'] = 0.75
        if parameters['predator'] == 'adapting':
            logger.info('Setting adaptation predator to begin (adaptbegin) at genrep*{0}'.format(parameters['adaptbegin']))
    if 'adaptmultiplier' not in parameters:
        parameters['adaptmultiplier'] = 3.0
        if parameters['predator'] == 'adapting':
            logger.info('Setting adaptation predator multiplier (adaptmultiplier) = {0}'.format(parameters['adaptmulitplier']))
    if 'demin' not in parameters:
        parameters['demin'] = 0.005
        logger.info('Setting cutoff convergence energy (demin) = {0}'.format(parameters['demin']))

    return parameters
