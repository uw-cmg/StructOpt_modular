import numpy as np

def random_three_vector(r):
    """Generates a random 3D unit vector (direction) with a
    uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """

    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = r * np.sin( theta) * np.cos( phi )
    y = r * np.sin( theta) * np.sin( phi )
    z = r * np.cos( theta )
    return (x,y,z)


def random_velocities(individual):
    velocities = np.zeros((len(individual), 3), dtype=np.float)
    for row in velocities:
        try:
            vx, vy, vz = random_three_vector(individual.pso_moves_parameters.update_particles.kwargs.vel_mag)
        except AttributeError:
            vx, vy, vz = random_three_vector(1)
        row[0], row[1], row[2] = vx, vy, vz
    individual.set_velocities(velocities)

