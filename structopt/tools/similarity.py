import numpy as np
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ase.cluster.octahedron import Octahedron
from ase.visualize import view
from structopt.cluster.individual.generators import fcc

def get_3d_psf(dimensions, r, HWHM):
    """Generates a psf array built from a gaussian function. The relevant 
    parameters specified in the parameters dictionary are below."""

    a, b, c = dimensions

    N_a = int(a * r)
    N_b = int(b * r)
    N_c = int(c * r)

    k_a_min = -0.5 * r
    k_a_max = (0.5 - 1.0/N_a) * r
    k_b_min = -0.5 * r
    k_b_max = (0.5 - 1.0/N_b) * r
    k_c_min = -0.5 * r
    k_c_max = (0.5 - 1.0/N_c) * r

    k_a = np.linspace(k_a_min, k_a_max, N_a)
    k_b = np.linspace(k_b_min, k_b_max, N_b)
    k_c = np.linspace(k_c_min, k_c_max, N_c)

    Mk_a, Mk_b, Mk_c = np.meshgrid(k_a, k_b, k_c)
    Mksq = Mk_a ** 2 + Mk_b ** 2 + Mk_c ** 2

    d_k =  (2 * np.log(2)) ** 0.5 / (HWHM * 2 * np.pi)

    psf = np.exp(-Mksq / (2 * d_k ** 2))

    return psf

def get_linear_convolution(dimensions, r, individual):
    """Calculate linear convoluted potential of an individual"""

    xmax, ymax, zmax = dimensions
    nx = int(xmax * r)
    ny = int(ymax * r)
    nz = int(zmax * r)
    dx = xmax/nx
    dy = ymax/ny
    dz = zmax/nz

    ax, ay, az = individual.get_positions().T

    # Assign atom to the bottom left of the grid point
    ix, iy, iz = np.floor(ax / dx), np.floor(ay / dy), np.floor(az / dz)

    # Apply periodic boundary conditions, considering all
    # corners of each pixel
    iax, ibx = np.fmod(ix, nx), np.fmod(ix + 1, nx)
    iay, iby = np.fmod(iy, ny), np.fmod(iy + 1, ny)
    iaz, ibz = np.fmod(iz, nz), np.fmod(iz + 1, nz)

    # Array of fraction of atoms at the left and bottom of the pixel
    fax = 1 - np.fmod(ax / dx, 1)
    fay = 1 - np.fmod(ay / dy, 1)
    faz = 1 - np.fmod(az / dz, 1)

    # Add potentials to grid. Split up the potential into
    # fractions on the pixel
    Zatom = np.ones(len(individual))
    V1 = fax * fay * faz * Zatom
    V2 = (1 - fax) * fay * faz * Zatom
    V3 = fax * (1 - fay) * faz * Zatom
    V4 = (1 - fax) * (1 - fay) * faz * Zatom
    V5 = fax * fay * (1 - faz) * Zatom
    V6 = (1 - fax) * fay * (1 - faz) * Zatom
    V7 = fax * (1 - fay) * (1 - faz) * Zatom
    V8 = (1 - fax) * (1 - fay) * (1 - faz) * Zatom

    V = np.zeros([nx, ny, nz])
    for j in range(len(individual)):
        V[int(iax[j]), int(iay[j]), int(iaz[j])] += V1[j]
        V[int(ibx[j]), int(iay[j]), int(iaz[j])] += V2[j]
        V[int(iax[j]), int(iby[j]), int(iaz[j])] += V3[j]
        V[int(ibx[j]), int(iby[j]), int(iaz[j])] += V4[j]
        V[int(iax[j]), int(iay[j]), int(iaz[j])] += V5[j]
        V[int(ibx[j]), int(iay[j]), int(iaz[j])] += V6[j]
        V[int(iax[j]), int(iby[j]), int(iaz[j])] += V7[j]
        V[int(ibx[j]), int(iby[j]), int(iaz[j])] += V8[j]

    return V
