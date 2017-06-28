import logging

import numpy as np

from snakebacon.bacon import run_baconmcmc

log = logging.getLogger(__name__)


class McmcProblem:
    def __init___(self, coredates, **kwargs):
        # TODO(brews): Write function to initialize McmcSetup
        pass

    def prior_dates(self):
        # TODO(brews): Write function for prior dates distributions
        pass

    def prior_sediment_rate(self):
        # TODO(brews): Write function for prior sediment rate distributions
        pass

    def prior_sediment_memory(self):
        # TODO(brews): Write function for prior sediment memory distributions
        pass


class McmcSetup(McmcProblem):
    def __init___(self, coredates, **kwargs):
        # TODO(brews): Write function to initialize McmcSetup
        pass

    def run(self):
        # TODO(brews): Write function to run McmcSetup.
        pass


class McmcResults(McmcProblem):
    def __init__(self, coredates, **kwargs):
        mcmcout = run_baconmcmc(core_labid=coredates.labid, core_age=coredates.age,
                                core_error=coredates.error, core_depth=coredates.depth, **kwargs)
        self.depth_segments = np.linspace(kwargs['depth_min'], kwargs['depth_max'], kwargs['k'])
        self.headage = mcmcout['theta']
        self.sediment_rate = mcmcout['x']
        self.sediment_memory = mcmcout['w']
        self.objective = mcmcout['objective']

    def burnin(self, n):
        """Remove the earliest n ensemble members from the MCMC output"""
        self.sediment_rate = self.sediment_rate[:, n:]
        self.headage = self.headage[n:]
        self.sediment_memory = self.sediment_memory[n:]
        self.objective = self.objective[n:]

    def n_members(self):
        """Get number of MCMC ensemble members in results"""
        return len(self.objective)

    def plot(self):
        # TODO(brews):
        pass
