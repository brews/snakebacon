from .bacon import run_baconmcmc


class Bacon:
    def runmcmc(*args, **kwargs):
        return run_baconmcmc(*args, **kwargs)

    def prior_dates(*args, **kwargs):
        """Get the prior distribution of calibrated radiocarbon dates"""
        pass

    def prior_sediment_rate(*args, **kwargs):
        """Get the prior distribution of sediment rates"""
        pass

    def prior_sediment_memory(*args, **kwargs):
        """Get the prior distribution of sediment memory"""
        pass
