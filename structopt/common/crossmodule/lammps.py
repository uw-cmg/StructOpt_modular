import os
import logging
import numpy as np
import shutil
import decimal as dec
from tempfile import mkdtemp, NamedTemporaryFile, mktemp as uns_mktemp
from re import compile as re_compile, IGNORECASE
from ase import Atoms, Atom

#from structopt.tools import root, single_core, parallel
from subprocess import Popen, PIPE, TimeoutExpired

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = '__end_of_ase_invoked_calculation__'

class LAMMPS(object):
    """Simplied calculator object for performing LAMMPS calculations
    through ase and StructOpt. Only returns the relaxed structure
    and the total energy. The primary difference between this and the
    ase version is that input files are written and then executed, while
    in the ase version the input is passed directly through the Popen
    constructor as stdin"""

    def __init__(self, parameters, calcdir=None):
        self.parameters = parameters
        self.cwd = os.getcwd()

        # read_log depends on that the first (three) thermo_style custom args
        # can be capitilized and matched against the log output. I.e.
        # don't use e.g. 'ke' or 'cpu' which are labeled KinEng and CPU.
        self._custom_thermo_args = ['step', 'temp', 'press', 'cpu',
                                    'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                                    'ke', 'pe', 'etotal',
                                    'vol', 'lx', 'ly', 'lz', 'atoms']
        self._custom_thermo_mark = ' '.join([x.capitalize() for x in
                                             self._custom_thermo_args[0:3]])

        # Match something which can be converted to a float
        f_re = r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))'
        n = len(self._custom_thermo_args)
        # Create a re matching exactly N white space separated floatish things
        self._custom_thermo_re = re_compile(r'^\s*' + r'\s+'.join([f_re]*n) + r'\s*$',
                                            flags=IGNORECASE)
        # thermo_content contains data "written by" thermo_style.
        # It is a list of dictionaries, each dict (one for each line
        # printed by thermo_style) contains a mapping between each
        # custom_thermo_args-argument and the corresponding
        # value as printed by lammps. thermo_content will be
        # re-populated by the read_log method.
        self.thermo_content = []
        self.pea = []

        # Initialize and set some default parameters here
        if calcdir:
            self.calcdir = calcdir
        elif 'calcdir' in self.parameters:
            self.calcdir = self.parameters['calcdir']
        else:
            self.calcdir = os.getcwd()

        self.parameters.setdefault('thermosteps', 0)
        self.parameters.setdefault('timeout', 60)
        self.parameters.setdefault('relax_box', False)

        return

    def calculate(self, atoms):
        self.atoms = atoms

        # Set periodic boundary conditions and get a prism-like
        # cell for running calculations
        cell = self.atoms.get_cell()
        self.prism = prism(cell)

        # Initialize and run the calculation
        self.setup_dir()
        os.chdir(self.tmp_dir)
        self.write_data()
        self.initialize()
        self.write_input()
        errors = self.run()

        if errors:
            self.process_error()

        # Read the thermodynamic and atom data
        self.read_log_file()
        self.read_trj_file()

        os.chdir(self.cwd)
        
        if self.parameters['keep_files'] == True:
            if not os.path.isdir(self.calcdir):
                os.makedirs(self.calcdir)
            for f in os.listdir(self.tmp_dir):
                f = os.path.join(self.tmp_dir, f)
                shutil.copy(f, self.calcdir)

        shutil.rmtree(self.tmp_dir)

        return

    def setup_dir(self):
        """This function sets up the temporary directory and copies
        the necessary files to that folder"""
        
        self.tmp_dir = mkdtemp(prefix='LAMMPS-')
        for param in self.parameters:
            if param.endswith('_file'):
                f = os.path.expandvars(self.parameters[param])
                shutil.copy(f, os.path.join(self.tmp_dir, os.path.basename(f)))

        self.trj_file = os.path.join(self.tmp_dir, 'trj.lammps')

        return

    def write_data(self):
        """Function for writing the atom positions in a seperate file"""

        atoms = self.atoms
        prism = self.prism

        self.data_file = os.path.join(self.tmp_dir, 'data.lammps')
        f = open(self.data_file, 'w')

        f.write('{} (written by ASE) \n\n'.format(f.name))
        atoms.wrap()
        atoms.center()
        symbols = atoms.get_chemical_symbols()
        n_atoms = len(symbols)
        f.write('{} \t atoms \n'.format(n_atoms))
        species = sorted(set(symbols))
        n_atom_types = len(species)
        f.write('{}  atom types\n'.format(n_atom_types))

        pbc = self.atoms.get_pbc()
        xhi, yhi, zhi, xy, xz, yz = prism.get_lammps_prism_str()
        xyzhis = [xhi, yhi, zhi]
        for index, axis in enumerate(['x','y','z']):
            if pbc[index]:    
                f.write('0.0 {}  {}lo {}hi\n'.format(xyzhis[index], axis, axis))
            else:
                xlo = min([ self.atoms.get_positions()[id][index] for id in range(len(self.atoms.get_positions())) ])
                xhi = max([ self.atoms.get_positions()[id][index] for id in range(len(self.atoms.get_positions())) ])
                f.write('{} {}  {}lo {}hi\n'.format(xlo, xhi, axis, axis))
        
        if prism.is_skewed():
            f.write('{} {} {}  xy xz yz\n'.format(xy, xz, yz))
        
        f.write('\n\n')

        f.write('Atoms \n\n')
        for i, r in enumerate(map(prism.pos_to_lammps_str, atoms.get_positions())):
            s = species.index(symbols[i]) + 1
            line = '{:>6} {:>3} {} {} {}\n'
            f.write(line.format(*(i+1, s)+tuple(r)))

        f.close()
        return

    def write_input(self):
        """Method which writes the LAMMPS in file"""

        self.input_file = os.path.join(self.tmp_dir, 'input.lammps')
        f = open(self.input_file, 'w')

        f.write('# (written by ASE)\n')

        # Write variables
        f.write('clear\n')
        f.write('variable dump_file string "{}"\n'.format(self.trj_file))
        f.write('variable data_file string "{}"\n'.format(self.data_file))

        # Write the atoms data
        parameters = self.parameters
        pbc = self.atoms.get_pbc()
        f.write('units metal \n')
        f.write('boundary {} {} {} \n'.format(*('sp'[x] for x in pbc)))
        f.write('read_data {}\n'.format(self.data_file))

        # Write interaction parameters
        f.write('\n### interactions \n')
        for param in ['pair_style', 'pair_coeff', 'mass']:
            if param in parameters:
                f.write('{} {}\n'.format(param, parameters[param]))

        # Write thermo parameters
        f.write('thermo_style custom {}\n'.format(' '.join(self._custom_thermo_args)))
        f.write('thermo_modify flush yes\n')
        f.write('thermo {}\n'.format(parameters['thermosteps']))

        # Relax the system
        f.write('\n### Relaxation \n')
        f.write('fix fix_nve all nve\n')
        if parameters['relax_box']:
            f.write('fix relax_box all box/relax iso 0.0 vmax 0.001\n')
        for param in ['min_style', 'min_modify', 'minimize']:
            if param in parameters:
                f.write('{} {}\n'.format(param, parameters[param]))
        f.write('compute pea all pe/atom\n')

        # Generate the thermodynamic and structural information
        dump_line = 'dump dump_all all custom 2 {} id type x y z c_pea\n'
        f.write(dump_line.format(self.trj_file))
        f.write('run 1\n')
        f.write('print {}'.format(CALCULATION_END_MARK))

        f.close()
        
        return

    def initialize(self):
        """The purpose of this function is to initialize the variables for 
        them being written to the lammps. This set of initializations
        depends on the atoms object, and hence cannot be done in __init__.py."""

        atoms = self.atoms

        # Initialize the potential parameters
        if 'pair_style' not in self.parameters:
            self.parameters['pair_style'] = 'lj/cut 10.0'
            self.parameters['pair_coeff'] = '* * 1 1'
            self.parameters['mass'] = '* 1.0'

        elif self.parameters['pair_style'] == 'eam':
            pot_file = os.path.expandvars(self.parameters['potential_file'])
            self.parameters['pair_coeff'] = '* * {}'.format(pot_file)

        elif self.parameters['pair_style'] == 'eam/alloy':
            elements = sorted(set(atoms.get_chemical_symbols()))
            pot_file = os.path.expandvars(self.parameters['potential_file'])
            pair_coeff = '* * {}'.format(pot_file)
            for element in elements:
                pair_coeff += ' {}'.format(element)
            self.parameters['pair_coeff'] = pair_coeff
        elif self.parameters['pair_style'] == 'eam/fs':
            elements = sorted(set(atoms.get_chemical_symbols()))
            pot_file = os.path.expandvars(self.parameters['potential_file'])
            pair_coeff = '* * {}'.format(pot_file)
            for element in elements:
                pair_coeff += ' {}'.format(element)
            self.parameters['pair_coeff'] = pair_coeff  
        elif 'lj/cut' in self.parameters['pair_style']:
            self.parameters['pair_coeff'] = '* * 1 1'
            self.parameters['mass'] = '* 1.0'
        else:
            s = self.parameters['pair_style']
            raise NotImplementedError('{} pair_style not yet implemented'.format(s))

        return

    def run(self):
        if 'LAMMPS_COMMAND' in os.environ:
            lammps_cmd_line = os.environ['LAMMPS_COMMAND']
        else:
            os.chdir(self.cwd)
            shutil.rmtree(self.tmp_dir)
            raise RuntimeError('Please set LAMMPS_COMMAND environment variable')

        input_file = open(self.input_file)
        p = Popen([lammps_cmd_line], stdin=input_file, stdout=PIPE, stderr=PIPE)
        try:
            output, error = p.communicate(timeout=self.parameters['timeout'])
        except TimeoutExpired:
            print("Timed out!")    
            return True

        self.output = output.decode('utf-8').split('\n')[:-1]

        # Check if the calculation completed without errors. If it does,
        # we need to save the files self.calcdir.
        if len(self.output) == 0 or CALCULATION_END_MARK not in self.output[-1]:
            return True

        return False

    def read_log_file(self):
        """Method which reads a LAMMPS output log file. This reads exclusively
        for the thermodynamic data."""

        if hasattr(self, 'output'):
            lines = self.output
        elif self.parameters['keep_files'] == True:            
            with open('{}/log.lammps'.format(self.calcdir)) as f:
                lines = f.readlines()
        else:
            raise RuntimeError('No log file detected. ' 
                               'Calculation not run or output not saved')

        thermo_content = []
        reading_thermo = False
        for line in lines:
            # get thermo output
            if line.startswith(self._custom_thermo_mark):
                reading_thermo = True
                continue

            thermo_step = self._custom_thermo_re.match(line)
            if reading_thermo and not thermo_step:
                reading_thermo = False
                continue            
            elif reading_thermo:
                # create a dictionary between each of the thermo_style args
                # and it's corresponding value
                thermo_content.append(dict(zip(self._custom_thermo_args,
                                               map(float, thermo_step.groups()))))

        self.thermo_content = thermo_content
        self.energy = thermo_content[-1]['pe']

        return

    def read_trj_file(self):
        """Method which reads the LAMMPS trj file. This is read primarily
        to get the atoms final relaxed structure"""

        if os.path.basename(os.getcwd()) == os.path.basename(self.tmp_dir):
            with open('trj.lammps') as f:
                lines = f.readlines()
        elif self.parameters['keep_files'] == True:
            with open('{}/log.lammps'.format(self.calcdir)) as f:
                lines = f.readlines()
        else:
            raise RuntimeError('No trajectory file detected. '
                               'Calculation not run or output not saved')

        # Get a list referencing atoms to lammps types
        atoms = self.atoms
        species = sorted(set(atoms.get_chemical_symbols()))

        for i, line in enumerate(lines):

            if 'ITEM: TIMESTEP' in line:
                lo = [] ; hi = [] ; tilt = []
                id = [] ; type = []
                positions = [] ; pea = [] #; velocities = [] ; forces = []

            if 'ITEM: NUMBER OF ATOMS' in line:                
                n_atoms = int(lines[i + 1].split()[0])

            if 'ITEM: BOX BOUNDS' in line:
                tilt_items = line.split()[3:]
                for j in range(3):
                    box_line = lines[i + j + 1]                    
                    fields = box_line.split()
                    lo.append(float(fields[0]))
                    hi.append(float(fields[1]))
                    if (len(fields) >= 3):
                        tilt.append(float(fields[2]))

            if 'ITEM: ATOMS' in line:
                atom_lines = [l.split() for l in lines[i+1:i+1+n_atoms]]
                ids, types, xs, ys, zs, peas = zip(*atom_lines)
                syms = [species[int(i) - 1] for i in types]
                pos = [None for i in range(len(ids))]
                peas = [None for i in range(len(ids))]
                for id, x, y, z, E in zip(ids, xs, ys, zs, peas):
                    pos[int(id) - 1] = [float(x), float(y), float(z)]
                    #peas[int(id) - 1] = float(E)

                # Update the positions of the atom
                self.atoms.set_positions(pos)

        # determine cell tilt (triclinic case!)
        if (len(tilt) >= 3):
            if (len(tilt_items) >= 3):
                xy = tilt[tilt_items.index('xy')]
                xz = tilt[tilt_items.index('xz')]
                yz = tilt[tilt_items.index('yz')]
            else:
                xy = tilt[0]
                xz = tilt[1]
                yz = tilt[2]
        else:
            xy = xz = yz = 0
        xhilo = (hi[0] - lo[0]) - xy - xz
        yhilo = (hi[1] - lo[1]) - yz
        zhilo = (hi[2] - lo[2])

        cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]]
        if all(atoms.get_pbc()):
            self.atoms.set_cell(cell)
                
        return
        
    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def update(self, atoms):
        if not hasattr(self, 'atoms') or self.atoms != atoms:
            self.calculate(atoms)

    def process_error(self):
        """This function is run immediately after detecting an error.
        We're in the temporary directory, so we have to copy the files
        back to the calculation directory, write and empty error file
        and raise an exception"""

        if not os.path.isdir(self.calcdir):
            os.makedirs(self.calcdir)
        for f in os.listdir(self.tmp_dir):
            f = os.path.join(self.tmp_dir, f)
            if os.path.isfile(f):
                shutil.copy(f, self.calcdir)
            elif os.path.isdir(f):
                shutil.move(f, self.calcdir)
            else:
                raise ValueError("The thing trying to be copied is not a file or directory")

        error_file = os.path.join(self.calcdir, 'error')
        open(error_file, 'a').close()
        os.chdir(self.cwd)

        raise RuntimeError('Error in LAMMPS calculation in {}'.format(self.calcdir))

class prism:
    def __init__(self, cell, pbc=(True,True,True), digits=10):
        """Create a lammps-style triclinic prism object from a cell

        The main purpose of the prism-object is to create suitable
        string representations of prism limits and atom positions
        within the prism.
        When creating the object, the digits parameter (default set to 10)
        specify the precission to use.
        lammps is picky about stuff being within semi-open intervals,
        e.g. for atom positions (when using create_atom in the in-file),
        x must be within [xlo, xhi).
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]
        alpha = np.arccos(np.dot(b, c)/(bn*cn))
        beta  = np.arccos(np.dot(a, c)/(an*cn))
        gamma = np.arccos(np.dot(a, b)/(an*bn))
        
        xhi = an
        xyp = np.cos(gamma)*bn
        yhi = np.sin(gamma)*bn
        xzp = np.cos(beta)*cn
        yzp = (bn*cn*np.cos(alpha) - xyp*xzp)/yhi
        zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)
    
        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
            int(np.floor(np.log10(max((xhi,yhi,zhi))))-digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        Apre = np.array(((xhi, 0,   0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), Apre)

        # Actual lammps cell may be different from what is used to create R
        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5*p
            n = (np.mod(x, p) - x)/p
            return [float(self.f2qdec(a)) for a in (vec + n*pvec)]

        Apre[1,:] = fold(Apre[1,:], Apre[0,:], 0)
        Apre[2,:] = fold(Apre[2,:], Apre[1,:], 1)
        Apre[2,:] = fold(Apre[2,:], Apre[0,:], 0)

        self.A = Apre
        self.Ainv = np.linalg.inv(self.A)

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def f2qs(self, f):
        return str(self.f2qdec(f))

    def f2s(self, f):
        return str(dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_HALF_EVEN))

    def dir2car(self, v):
        "Direct to cartesian coordinates"
        return np.dot(v, self.A)

    def r2dir(self, v):
        "Cartesian to direct coordinates"
        return np.dot(v, self.Ainv)

    def fold_to_str(self,v):
        "Fold a position into the lammps cell (semi open), return a tuple of str"
        # Two-stage fold, first into box, then into semi-open interval
        # (within the given precission).
        d = [x % (1-self.dir_prec) for x in
             map(dec.Decimal, map(repr, np.mod(self.car2dir(v) + self.eps, 1.0)))]
        return tuple([self.f2qs(x) for x in
                      self.dir2car(list(map(float, d)))])
        
    def get_lammps_prism(self):
        A = self.A
        return (A[0,0], A[1,1], A[2,2], A[1,0], A[2,0], A[2,1])

    def get_lammps_prism_str(self):
        "Return a tuple of strings"
        p = self.get_lammps_prism()
        return tuple([self.f2s(x) for x in p])

    def pos_to_lammps_str(self, position):
        "Rotate an ase-cell position to the lammps cell orientation, return tuple of strs"
        return tuple([self.f2s(x) for x in np.dot(position, self.R)])

    def pos_to_lammps_fold_str(self, position):
        "Rotate and fold an ase-cell position into the lammps cell, return tuple of strs"
        return self.fold_to_str(np.dot(position, self.R))

    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)
