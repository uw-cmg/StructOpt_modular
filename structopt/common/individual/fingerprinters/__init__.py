import functools

import structopt
from . import dirac
from structopt.tools import root, single_core, parallel


class Fingerprinters(object):
    """ """

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.fingerprinters
        self.fingerprinters = [getattr(self, name) for name in self.parameters.options]
        self.selected_fingerprinter = None

    @single_core
    def select_fingerprinter(self):
        self.selected_fingerprinter = random.choice(self.fingerprinters)

    @single_core
    def fingerprint(self, individual):
        return self.selected_fingerprinter(individual)

    @single_core
    def post_processing(self):
        pass

    @staticmethod
    @functools.wraps(dirac)
    def dirac(individual):
        return dirac(individual)

