import logging

from snakebacon.bacon import run_baconmcmc

log = logging.getLogger(__name__)


class McmcResults:
    def __init__(self, coredates, **kwargs):
        # mcmcout = run_baconmcmc(core_labid = coredates.labid, core_age = coredates.age,
        #                                   core_error = coredates.error, core_depth = coredates.depth,
        #                                   depth_min = 1.5, depth_max = 99.5, cc=[1],
        #                                   cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4 = 'ConstCal',
        #                                   d_r = [0], d_std = [0], t_a=[3], t_b=[4], k= 20,
        #                                   minyr=-1000, maxyr = 1e6, th01=4147, th02=4145,
        #                                   acc_mean = 20, acc_shape = 1.5, mem_strength = 4, mem_mean = 0.7)
        mcmcout = run_baconmcmc(core_labid=coredates.labid, core_age=coredates.age,
                                core_error=coredates.error, core_depth=coredates.depth, **kwargs)
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
        pass
