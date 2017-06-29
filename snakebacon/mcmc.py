import logging

import numpy as np

from snakebacon.records import DateRecord
from snakebacon.bacon import run_baconmcmc

log = logging.getLogger(__name__)


class McmcProblem:
    def __init__(self, coredates, **kwargs):
        self.coredates = DateRecord(coredates)
        self.mcmc_kwargs = kwargs

    def prior_dates(self):
        """Get the prior distribution of radiocarbon dates"""
        # TODO(brews): Write function for prior dates distributions
        pass

    def prior_sediment_rate(self):
        """Get the prior distribution of sediment rates"""
        # TODO(brews): Write function for prior sediment rate distributions
        pass

    def prior_sediment_memory(self):
        """Get the prior distribution of sediment memory"""
        # TODO(brews): Write function for prior sediment memory distributions
        pass


class McmcSetup(McmcProblem):
    def __init___(self, coredates, **kwargs):
        super().__init__(coredates, **kwargs)

    def run(self):
        return McmcResults(self.coredates, **self.mcmc_kwargs)


class McmcResults(McmcProblem):
    def __init__(self, coredates, **kwargs):
        super().__init__(coredates, **kwargs)
        mcmcout = run_baconmcmc(core_labid=self.coredates.labid, core_age=self.coredates.age,
                                core_error=self.coredates.error, core_depth=self.coredates.depth, **self.mcmc_kwargs)
        self.depth_segments = np.linspace(self.mcmc_kwargs['depth_min'], self.mcmc_kwargs['depth_max'], self.mcmc_kwargs['k'])
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
        # TODO(brews): Write function to plot raw Mcmc results.
        pass
