import logging

import numpy as np

from snakebacon.records import DateRecord
from snakebacon.bacon import run_baconmcmc


log = logging.getLogger(__name__)


class McmcSetup:
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

    def validate(self):
        """Validate configs for mcmc run"""
        # TODO(brews): Write mcmc config validation
        pass

    def run(self):
        self.validate()
        return McmcResults(self)


class McmcResults:
    def __init__(self, setup):
        mcmcout = run_baconmcmc(core_labid=setup.coredates.labid, core_age=setup.coredates.age,
                                core_error=setup.coredates.error, core_depth=setup.coredates.depth,
                                **setup.mcmc_kwargs)
        self.depth_segments = np.linspace(setup.mcmc_kwargs['depth_min'], setup.mcmc_kwargs['depth_max'],
                                          setup.mcmc_kwargs['k'])
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
