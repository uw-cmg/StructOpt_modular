'''This contains a StructOpt calculator for submitting jobs to queue, tracking their progress, and reading their output. Currently only works with genetic.py'''

import os

class StructOpt(object):

    def __init__(self, calcdir=None, optimizer=None, parameters=None):
        if calcdir == None:
            self.
