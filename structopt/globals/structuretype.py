import logging
import random

import structopt


class StructureType(object):
	"""An abstract base class for a structure type."""

	def __init__(self, cls_name):
		self.logger = logging.getLogger('default')
		for k, v in structopt.parameters.items():
            setattr(self, k, v)

        self.relaxation_modules = []
        for relaxation in self.relaxations:
            mod = import_module('structopt.{cls_name}.relaxations.{module_name}'.format(cls_name=cls_name, module_name=relaxation))  # Import the module package
            cls = getattr(mod, '{cls_name}_eval'.format(cls_name=relaxation))  # Get's the class from the module package
            self.relaxation_modules.append(cls())

        self.fitness_modules = []
        for fitness in self.modules:
            mod = import_module('structopt.{cls_name}.fitnesses.{module_name}'.format(cls_name=cls_name, module_name=fitness))  # Import the module package
            cls = getattr(mod, '{cls_name}_eval'.format(cls_name=fitness))  # Get's the class from the module package
            self.fitness_modules.append(cls())

    	# Setup the output files

    	# Set starting convergence and generations

    	# Prep output monitoring

    	# Initialize random number seed
    	random.seed(self.seed)

    	# Write the input parameters to the output file
    	self.logger.debug('Writing the input parameters to output file')
        structopt.fileio.parameters.write(self)

        # Allow structopt to see me
        structopt.parameters.globals.structuretype = self
