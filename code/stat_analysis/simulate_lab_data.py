""" 
Description
-----------

Simulating experiment to test model

Want to check that the latent zoospore model we are using can approximately recover 
known parameter values. The below code sets up a simulation of the mesocosm 
experiment with known parameters and then fits a Stan model to estimate the 
parameters.  We also test whether the latent zoospore model can 
correctly distinguish between two transmission functions.

Case 1: Ideal Case
-------------------
We have frog and tadpole measurements every day for 32 days

Verdict: Model can recover the transmission parameters.

Case 2: Longer between swab times
---------------------------------
We have frog and tadpole measurements every 4 days.

Verdict: Model can recover the transmission parameters. 

Case 3: Unequal time between swabs
----------------------------------

Verdict: Model can recover the parameters even with unequal swab times.
There might be a slight bias in the parameters, particularly towards
underestimating the important of transmission.

Case 4: Compare DIC values between true and false models
-----------------------------------

Calculated DIC from a density-dependent model (true) and a 
frequency-dependent model (false) and the DIC for the true model much lower.

This indicates that this approach can identify the parameters and distinguish
between models.


"""

from __future__ import division
import numpy as np
import pandas as pd
import pystan
import scipy.stats as stats


class Frog(object):
    """ 
    A frog in a tank.  Just as is assumed in the host-parasite IPM,
    a frog can get infected with Bd with some initial load, that load can grow,
    and the frog can lose that load.  All of these are probabilistic events 
    with the same functional forms as given in the frog-Bd IPM described in the
    manuscript.

    Attributes
    ----------
    : see below
    """

    def __init__(self, load, growth_params):

        self.growth_params = growth_params
        self.load = load

    def is_infected(self):
        return(self.load > 0)

    def get_loss(self):
        # Loss function of IPM

        coef = np.exp(self.growth_params['loss_int'] + 
                      self.growth_params['loss_slope']*np.log(self.load))
        loss_prob = coef / (1 + coef)
        loss = np.random.binomial(1, loss_prob, 1)
        return(loss)

    def get_bd_growth(self):
        # Bd growth function of IPM

        mu = self.growth_params['growth_int'] + \
                    self.growth_params['growth_slope'] * np.log(self.load)
        new_load = np.random.normal(loc=mu, scale=growth_params['growth_sd'], 
                                            size=1)
        return(new_load)

    def get_transmission(self, num_inf, zpool):
        # Transmission function of IPM

        inner = self.growth_params['inf_int']*zpool + \
                self.growth_params['inf_slope']*num_inf
        prob_inf = 1 - np.exp(-(inner))
        inf = np.random.binomial(1, prob_inf, size=1)

        return(inf)

    def get_intial_load(self):
        # Initial load of infection function of IPM

        init_load = np.random.normal(loc=self.growth_params['init_inf_mu'],
                                  scale=self.growth_params['init_inf_sd'],
                                  size=1)
        return(init_load)

    def update_bd_load(self, num_inf, zpool):
        # Update bd load on a frog for a single time step

        if self.is_infected(): # Change load if infected

            # Lose infection
            loss = self.get_loss()
            new_load = self.get_bd_growth()

            self.load = float((1 - loss) * np.exp(new_load))

        else: # Determine if the frog becomes infected

            inf = self.get_transmission(num_inf, zpool)
            init_load = self.get_intial_load()
            
            self.load = float(inf * np.exp(init_load))

        return(self.load)


class Tank(object):
    """ A tank holds some number of frogs """

    def __init__(self, num_frogs, zpool_init, growth_params):

        self.num_frogs = num_frogs

        self.frogs = [Frog(0, growth_params) for i in range(num_frogs)]
        self.zpool = zpool_init
        self.growth_params = growth_params

    def update_zpool(self, frog_loads):
        """ 
        Update the zpool pool in the tank

        Parameters
        -----------
        float : frog_loads
            The summed Bd load on frogs in the tank at a time point

        """

        mu = self.zpool*np.exp(-self.growth_params['zdeath']) + frog_loads + 2000

        zpool_next = np.random.normal(loc=np.log(mu) - 0.5, 
                                      scale=self.growth_params['zsd'],
                                      size=1)[0]
        return(zpool_next)


    def update_tank(self):
        """ Update the frogs and zoospore pool in a tank"""

        # Update zoospore pool
        prev_zpool = self.zpool
        frog_loads = np.sum([frog.load for frog in self.frogs])
        self.zpool = np.exp(self.update_zpool(frog_loads))

        # Update frog loads
        num_inf = np.sum([frog.is_infected() for frog in self.frogs])
        trajs = [frog.update_bd_load(num_inf, np.log(prev_zpool + 1)) for frog in self.frogs]
        return((self.zpool, trajs))

def extend_array(array, max_row):
    """
    Extend array to have at least max_row rows
    """

    shape = array.shape
    diff = max_row - shape[0]

    if diff != 0:
        new_array = np.row_stack([array, np.full((diff, shape[1]), -1, dtype=np.int)])
    else:
        new_array = array

    return(new_array)


if __name__ == "__main__":

    # IPM parameter for simulation
    growth_params = {'growth_int': 1.577, "growth_slope": 0.799, "growth_sd": np.sqrt(5.92),
                     'loss_int': -1.35, 'loss_slope': -0.472,
                     'inf_int': 0.00529, 'inf_slope': 0.0752,
                     'init_inf_mu': 1.135, 'init_inf_sd': np.sqrt(0.59),
                     'zdeath': 0.3, 'zsd': 1}
    time_arrays = pd.read_csv("time_arrays.csv") # Load in empirical unequal swab times

    # Experimental parameters               
    frogs_per_tank = np.repeat([16, 8, 4, 1], 4)
    T = len(frogs_per_tank)
    N = 32 # Number of time points
    M = np.max(frogs_per_tank)

    time_step = 4 # Time between swabs, i.e. 4 would be swabbing every 4 days
    steps = np.arange(N + 1, step=time_step) # Time step indices
    unequal_swab = False # If this is True time_step must be 1

    if unequal_swab:
        time_step = 1

        
    prob_inf = []
    frog_zes = []
    density = []
    prev = []
    tad_zes = []
    time = []

    for num_frogs in frogs_per_tank:

        tank = Tank(num_frogs, 0, growth_params)

        # Run experiment for N days and extract data
        trajs = [tank.update_tank() for i in range(N)]
        ztraj, frog_trajs = zip(*trajs)

        # Format and add initial values add on initial values
        frog_dat = pd.DataFrame(np.array(frog_trajs)).T
        frog_dat = pd.DataFrame(np.hstack([np.repeat(0, num_frogs).reshape(num_frogs, 1), 
                              frog_dat]))


        prob_inf.append(np.array(extend_array((frog_dat > 0).astype(np.int), 
                                                            M))[:, steps])

        frog_zes.append(np.array(frog_dat.sum(axis=0))[steps])
        density.append(np.array((frog_dat > 0).sum(axis=0))[steps])
        prev.append((np.array((frog_dat > 0).sum(axis=0)) / frog_dat.shape[0])[steps])
        tad_zes.append(np.array(np.exp(np.random.normal(loc=7.49, 
                                            scale=1.65, size=len(steps)))))
        time.append(np.repeat(time_step, len(steps)))


    # If you have unequal swabbing times, drop data points based on the observed
    # time differences in the mesocosm experiment data
    if unequal_swab:
        cum_times = time_arrays.cumsum(axis=1)

        prob_inf = [prob_inf[i][:, row] for i, row in enumerate(np.array(cum_times))]
        frog_zes = [frog_zes[i][row] for i, row in enumerate(np.array(cum_times))]
        density = [density[i][row] for i, row in enumerate(np.array(cum_times))]
        prev = [prev[i][row] for i, row in enumerate(np.array(cum_times))]
        tad_zes = [tad_zes[i][row] for i, row in enumerate(np.array(cum_times))]
        time = [t for t in np.array(time_arrays)]
        steps = time[0] # just for the right length

    # Host many points are contributing to the likelihood?
    counter = 0
    for t in range(T):
        for i in range(len(steps)):
            if i != 0:
                for j in range(M):
                    if j <= frogs_per_tank[t]:
                        if prob_inf[t][j, i - 1] == 0:
                            counter = counter + 1


    stan_data = {'N': len(steps), 
                 'M': M, 
                 'T': T,
                 'num_frogs': frogs_per_tank,
                 'frog_zes': frog_zes,
                 'tad_zes': tad_zes,
                 'prob_inf': prob_inf,
                 'density': density,
                 'time': time,
                 'init_tads': 2000,
                 'D': counter,
                 'sigma': 1}

    mod_fit = pystan.stan(file='stan_files/pomp_zpool_with_decay.stan', 
                data=stan_data, warmup=2000, iter=6000, chains=3,
                n_jobs=3)

    stan_data_fd = {'N': len(steps), 
                 'M': M, 
                 'T': T,
                 'num_frogs': frogs_per_tank,
                 'frog_zes': frog_zes,
                 'tad_zes': tad_zes,
                 'prob_inf': prob_inf,
                 'density': prev,
                 'time': time,
                 'init_tads': 2000,
                 'D': counter,
                 "sigma": 1}

    mod_fit_fd = pystan.stan(file='stan_files/pomp_zpool_with_decay.stan', 
                data=stan_data_fd, warmup=2000, iter=6000, chains=3,
                n_jobs=3)


    def dic(mod):
        # DIC estimate for stan model

        dev_full = mod.extract()['dev']
        dev_bar = np.mean(dev_full)
        pd_alt = 2*np.var(dev_full, ddof=1)
        return(dev_bar + pd_alt)


    # Density-dependent model is better under DIC
    dics_vals = [dic(mod) for mod in (mod_fit, mod_fit_fd)]

    # Extract parameter estimates
    bzpool_ests = stats.scoreatpercentile(mod_fit.extract()['b_zpool'], 
                            (2.5, 50, 97.5))
    bdens_ests = stats.scoreatpercentile(mod_fit.extract()['b_density'], 
                            (2.5, 50, 97.5))



