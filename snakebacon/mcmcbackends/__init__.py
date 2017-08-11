from .bacon import run_baconmcmc


class Bacon:
    def runmcmc(*args, **kwargs):
        return run_baconmcmc(*args, **kwargs
                             )