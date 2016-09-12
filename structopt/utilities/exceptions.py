class StructOptUnknownState(Exception):
    pass

class StructOptRunning(Exception):
    pass

class StructOptQueued(Exception):
    pass

class StructOptSubmitted(Exception):
    def __init__(self, jobdir):
        self.jobdir = jobdir
    def __str__(self):
        return repr(self.jobdir)
