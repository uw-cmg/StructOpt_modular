import functools.wraps

import structopt.parameters
from . import dirac


class Fingerprinters(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.fingerprinters
        self.fingerprinters = [getattr(self, name) for name in self.parameters.fingerprinters.options]
        self.selected_fingerprinter = None

    def select_fingerprinter(self):
        self.selected_fingerprinter = random.choice(self.fingerprinters)

    def fingerprint(self, individual):
        return self.selected_fingerprinter(individual)

    @staticmethod
    @functools.wraps(dirac)
    def dirac(individual):
        return dirac(individual)

