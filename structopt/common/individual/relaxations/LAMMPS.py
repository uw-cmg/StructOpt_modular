class LAMMPS(object):
    def get_command(self, individual):
        return 'who'

    def get_energy(self, individual):
        return 0.0

    def relax(self, individual):
        #subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)
        return individual

