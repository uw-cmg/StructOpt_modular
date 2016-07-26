'''This contains a StructOpt calculator for submitting jobs to queue, tracking their progress, and reading their output. Currently only works with genetic.py'''

import os

class StructOpt(object):

    def __init__(self, calcdir=None, optimizer=None, parameters=None):

        if calcdir is None:
            self.calcdir = os.getcwd()
        else:
            self.calcdir = os.path.expanduser(calcdir)
        self.cwd = os.getcwd()
        self.optimizer = optimizer
        self.parameters = parameters

        self.system_name = os.path.basename(self.calcdir)

    def __enter__(self):
        '''On enter, make sure directory exists. Create it if necessary
        and change into the directory. Then return the calculator.'''

        # Make directory if it doesn't already exist
        if not os.path.isdir(self.calcdir):
            os.makedirs(self.calcdir)

        # Now change into the new working directory
        os.chdir(self.calcdir)
        self.initialize()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''On exit, change back to the original directory'''

        os.chdir(self.cwd)

        return

    def initialize(self):
        '''Gets the status of the calculation. Meant to be run when within
        the calculation directory'''

        pass
