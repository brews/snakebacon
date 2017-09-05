import logging

import numpy as np

import snakebacon.mcmcbackends
from snakebacon.records import ChronRecord


log = logging.getLogger(__name__)


class McmcSetup:
    def __init__(self, coredates, mcmcbackend=snakebacon.mcmcbackends.Bacon, **kwargs):
        self.coredates = ChronRecord(coredates)
        self.mcmc_kws = kwargs
        self.mcmcbackend = mcmcbackend

    def prior_dates(self):
        """Get the prior distribution of radiocarbon dates"""
        return self.mcmcbackend.prior_dates(self.coredates, **self.mcmc_kws)

    def prior_sediment_rate(self):
        """Get the prior distribution of sediment rates"""
        return self.mcmcbackend.prior_sediment_rate(self.coredates, **self.mcmc_kws)

    def prior_sediment_memory(self):
        """Get the prior distribution of sediment memory"""
        return self.mcmcbackend.prior_sediment_memory(self.coredates, **self.mcmc_kws)

    def validate(self):
        """Validate configs for mcmc run"""
        # TODO(brews): Write mcmc config validation
        pass

    def run(self):
        self.validate()
        return McmcResults(self)


class McmcResults:
    def __init__(self, setup):
        mcmcout = setup.mcmcbackend.runmcmc(core_labid=setup.coredates.labid,
                                            core_age=setup.coredates.age,
                                            core_error=setup.coredates.error,
                                            core_depth=setup.coredates.depth,
                                            **setup.mcmc_kws)
        self.depth_segments = np.linspace(setup.mcmc_kws['depth_min'], setup.mcmc_kws['depth_max'],
                                          setup.mcmc_kws['k'])
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
