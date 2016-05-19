import ase.io


def read_xyz(filename, index=None, format=None, **kwargs):
    return ase.io.read(filename, index, format, **kwargs)

