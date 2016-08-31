import os
import math
import logging
import numpy as np
from scipy.signal import fftconvolve

from ase.io import read

from structopt.tools import root, single_core, parallel
from structopt.tools.dictionaryobject import DictionaryObject

class STEM(object):
    """Calculates the chi^2 difference between a simulated and experimental image.
    In order to calculate a z-contrast image and chi^2 function the following
    parameters are required

    Parameters
    ----------
    HWHM : float
        The half-width half-maximum of the gaussian function used for the 
        point spread function.
    dimensions : list
        The x and y dimensions of STEM image.
    resolution : float
        The pixels per angstrom resolution
    """

    @single_core
    def __init__(self, parameters={}):
        self.parameters = parameters
        self.parameters.setdefault('zed', 1)
        self.psf = None
        self.target = None
        self.phantom = True

        # Set default normalization to E = E/nprotons
        if "normalize" not in self.parameters:
            self.parameters.setdefault("normalize", {})
        self.parameters['normalize'].setdefault('nprotons', True)
        self.parameters['normalize'].setdefault('SSE', False)

        # If running within StructOpt, create directory for saving files
        # and faster loading of PSF and target data
        if hasattr(logging, 'parameters'):
            self.path = os.path.join(logging.parameters.path, 'fitness/STEM')
            os.makedirs(self.path, exist_ok=True)
        else:
            self.path = None

    def fitness(self, individual):
        """Calculates the fitness of an individual with respect to a target
        image. Normalize this fitness by the number of atoms."""

        if self.psf is None:
            self.generate_psf()
        if self.target is None:
            self.generate_target()

        image = self.get_image(individual)
        image, x_shift, y_shift = self.cross_correlate(image)

        chi = image - self.target
        chi = self.normalize(chi, individual)

        return chi

    def cross_correlate(self, image):
        convolution = fftconvolve(image, self.target[::-1, ::-1], mode='full')
        y_max, x_max = np.unravel_index(np.argmax(convolution), convolution.shape)
        x_shift = x_max - image.shape[1] + 1
        y_shift = y_max - image.shape[0] + 1
        image = np.roll(image, x_shift, axis=1)
        image = np.roll(image, y_shift, axis=0)

        # Return the x_shift and y_shift in real space
        #x_shift /= self.parameters['resolution']
        #y_shift /= self.parameters['resolution']

        return image, x_shift, y_shift

    def normalize(self, chi, individual):
        if 'normalize' not in self.parameters:
            return chi

        norms = self.parameters['normalize']

        if 'SSE' in norms and norms['SSE']:
            chi = np.sum(np.square(chi)) ** 0.5
            if 'nprotons' in norms and norms['nprotons']:
                chi /= sum(individual.get_atomic_numbers()) ** 0.5
        else:
            chi = np.sum(np.absolute(chi))
            if 'nprotons' in norms and norms['nprotons']:
                chi /= sum(individual.get_atomic_numbers())

        return chi

    def get_Z_diff(self, individual):
        """Returns the difference between the individual and target's
        estimated atomic number count"""

        if self.target is None:
            self.generate_target()

        image = self.get_image(individual)
        zed = self.parameters['zed']
        image_Z_tot = np.sum(image ** (1/zed))
        if self.phantom == True:
            target_Z_tot = np.sum(self.target ** (1/zed))
        else:
            return NotImplementedError("Real STEM image not implemented yet")
        
        Z_diff = image_Z_tot - target_Z_tot

        return Z_diff

    def get_linear_convolution(self, individual):
        """Calculate linear convoluted potential of an individual"""

        r = self.parameters['resolution']
        zed = self.parameters['zed']
        xmax, ymax = self.parameters['dimensions']
        if type(xmax) is float:
            nx = int(xmax * r)
            ny = int(ymax * r)
            dx = xmax/nx
            dy = ymax/ny
        elif type(xmax) is int:
            nx = xmax
            ny = ymax
            xmax = nx / r
            ymax = ny / r
            dx = xmax / nx
            dy = ymax / ny

        ax, ay, az = individual.get_positions().T

        # Assign atom to the bottom left of the grid point
        ix, iy = np.floor(ax / dx), np.floor(ay / dy)
        
        # Apply periodic boundary conditions, considering all
        # corners of each pixel
        iax, ibx = np.fmod(ix, nx), np.fmod(ix + 1, nx)
        iay, iby = np.fmod(iy, ny), np.fmod(iy + 1, ny)

        # Array of fraction of atoms at the left and bottom of the pixel
        fax = 1 - np.fmod(ax / dx, 1)
        fay = 1 - np.fmod(ay / dy, 1)

        # Add potentials to grid. Split up the potential into
        # fractions on the pixel
        Zatom = individual.get_atomic_numbers()
        V1 = fax * fay * Zatom ** zed
        V2 = (1 - fax) * fay * Zatom ** zed
        V3 = fax * (1 - fay) * Zatom ** zed
        V4 = (1 - fax) * (1 - fay) * Zatom ** zed

        V = np.zeros([nx, ny])
        for j in range(len(individual)):
            V[iax[j],iay[j]] += V1[j]
            V[ibx[j],iay[j]] += V2[j]
            V[iax[j],iby[j]] += V3[j]
            V[ibx[j],iby[j]] += V4[j]

        return V

    def generate_psf(self):
        """Generates a psf array built from a gaussian function. The relevant 
        parameters specified in the parameters dictionary are below."""

        # We do not want to generate the psf if it has already been saved
        if (self.path is not None
            and os.path.isfile(os.path.join(self.path, 'psf.npy'))):
            self.psf = np.load(os.path.join(self.path, 'psf.npy'))
            print('psf loaded')
            return

        HWHM = self.parameters['HWHM']
        r = self.parameters['resolution']
        a, b = self.parameters['dimensions']

        if type(a) is float:
            N_a = int(a * r)
            N_b = int(b * r)
        else:
            N_a = a
            N_b = b

        k_a_min = -0.5 * r
        k_a_max = (0.5 - 1.0/N_a) * r
        k_b_min = -0.5 * r
        k_b_max = (0.5 - 1.0/N_b) * r

        k_a = np.linspace(k_a_min, k_a_max, N_a)
        k_b = np.linspace(k_b_min, k_b_max, N_b)

        Mk_a, Mk_b = np.meshgrid(k_a, k_b)
        Mksq = Mk_a ** 2 + Mk_b ** 2

        d_k =  (2 * np.log(2)) ** 0.5 / (HWHM * 2 * np.pi)

        psf = np.exp(-Mksq / (2 * d_k ** 2))

        self.psf = psf

        # Try saving the PSF for future calculations
        if self.path is not None:
            np.save(os.path.join(self.path, 'psf'), self.psf)

        return

    def generate_target(self):
        """Generates the target STEM image from the parameters. The bulk of
        this code generates the target from an atoms object. """

        # Load the target from a file
        if (self.path is not None
            and os.path.isfile(os.path.join(self.path, 'target.npy'))):
            self.target = np.load(os.path.join(self.path, 'target.npy'))
            print('target loaded')
            return

        if self.psf is None:
            self.generate_psf()

        if not self.parameters['target'].endswith('.xyz'):
            self.read_image(self.parameters['target'])
            self.process_image()
            return

        atoms = read(self.parameters['target'])
        self.target = self.get_image(atoms)
        self.phantom = True

        # Try saving the target for future calculations
        if self.path is not None:
            np.save(os.path.join(self.path, 'target'), self.target)

        return


    def get_image(self, individual):
        """Calculates the z-contrasted STEM image of an individual"""

        if self.psf is None:
            self.generate_psf()

        psf = self.psf
        V = self.get_linear_convolution(individual)

        ft_psf = np.fft.fftshift(psf)
        ft_V = np.fft.fft2(V).T

        zcon_image = np.fft.ifft2(ft_psf * ft_V, axes=(0, 1)).real

        return zcon_image


