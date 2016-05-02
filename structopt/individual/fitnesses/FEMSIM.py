import structopt.parameters


class FEMSIM(object):
    def __init__(self):
        self.parameters = structopt.paramters.fitnesses.FEMSIM

        femsimfiles = '{filename}-rank0/FEMSIMFiles'.format(filename=Optimizer.filename)
        if not os.path.exists(femsimfiles):
            os.mkdir(femsimfiles)

        self.vk = np.multiply(self.paramters.thickness_scaling_factor, self.vk)  # Multiply the experimental data by the thickness scaling factor

    def get_command(self):

    def get_chisq(self):

