import numpy as np
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ase.cluster.octahedron import Octahedron
from ase.visualize import view

from structopt.cluster.individual.generators import fcc
from structopt.tools.analysis import get_avg_radii

def get_chi2(atoms1, atoms2, cutoff=0.8, r=2.0, HWHM=0.4):
    """Calculates the chi2, which is the difference in positions
    between atoms1 and atoms2"""

    # Calculate the offset to apply to atoms1 to maximize matching
    # to atoms2
    offset = get_offset(atoms1, atoms2, r=r, HWHM=HWHM)
    atoms1.translate(offset)

    # Calculate the location difference of each atom1 atom with atom2 atom
    cutoff *= get_avg_radii(atoms1) * 2
    pos1 = np.expand_dims(atoms1.get_positions(), 0)
    pos2 = np.expand_dims(atoms2.get_positions(), 0)
    dists = np.linalg.norm(pos1 - np.transpose(pos2, (1, 0, 2)), axis=2)

    nn_dists = np.min(dists, axis=1)
    x_fp = np.sum(nn_dists > cutoff) / len(atoms2)
    x_fn = np.sum(np.min(dists, axis=0) > cutoff) / len(atoms2)
    chi2 = nn_dists[nn_dists <= cutoff]

    return x_fp, x_fn, chi2


def get_offset(atoms1, atoms2, r=5.0, HWHM=0.4):
    """Gets the offset to apply to atoms1 to have its positions match atoms2"""

    cell1 = atoms1.get_cell()
    cell2 = atoms2.get_cell()

    # Make sure the cell is a box and they are the same
    assert (cell1.diagonal() * np.eye(3) == cell1).all()
    assert (cell2.diagonal() * np.eye(3) == cell2).all()
    assert (cell1 == cell2).all()

    # Load the 3d point spread function (psf)
    psf = get_3d_psf(cell1.diagonal(), r, HWHM)
    ft_psf = np.fft.fftshift(psf)

    # Get the gridded locations
    V1 = get_gridded_locations(cell1.diagonal(), r, atoms1)
    ft_V1 = np.fft.fftn(V1).T
    image1 = np.fft.ifftn(ft_psf * ft_V1)

    V2 = get_gridded_locations(cell2.diagonal(), r, atoms2)
    ft_V2 = np.fft.fftn(V2).T
    image2 = np.fft.ifftn(ft_psf * ft_V2)

    # Use cross-correlation to calculate the ideal offset
    correlation = fftconvolve(image2, image1[::-1, ::-1, ::-1], mode='full')
    z_max, y_max, x_max = np.unravel_index(np.argmax(correlation), correlation.shape)
    x_shift = (x_max - image1.shape[2] + 1) / r
    y_shift = (y_max - image1.shape[1] + 1) / r
    z_shift = (z_max - image1.shape[0] + 1) / r

    return (x_shift, y_shift, z_shift)

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

def get_gridded_locations(dimensions, r, individual):
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
