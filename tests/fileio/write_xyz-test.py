
import sys; sys.path.insert(1, '../')
from structopt.fileio import write_xyz

from ase.lattice import bulk
from ase.visualize import view

atoms = bulk('Pt')
write_xyz('fileio/write_xyz-test.xyz', atoms, append=False)

with open('fileio/write_xyz-test.xyz') as f:
    print(f.read())
