class LAMMPS(object):
    def get_command(self, individual):
        return 'who'

    def get_energy(self, individual):
        return 0.0

    def relax(self, individual):
        #print("Running: {}".format(command))
        #subprocess.call(command)
        return individual

