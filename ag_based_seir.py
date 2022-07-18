import time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_backend
from acquaintance_network import acquaintance_network
from daily_contact_network import transmission_simulate


class AgentBasedSEIR:

    def __init__(self,
                 *,
                 # base params
                 nActors=20000,
                 percentage_infection=0.05,
                 househouldSize=3.25,
                 household_weight=3.35,
                 household_weight_outbreak=4.676,
                 geo_weight=0.5,
                 random_weight=0.5,
                 random_connect_num=1,
                 age_distribution_prob=None,
                 family_size=None,
                 family_age_distribution=None,
                 import_cases=None,
                 age_suscep=None,
                 remove_period=None,
                 outbreak_period=None,
                 window_exposed_ratio=None,
                 window_incubation=None,
                 window_post_infectious=None,
                 # scenario params
                 base_contact_matrix=None,
                 outbreak_contact_matrix=None,
                 min_remove=None,
                 household_quar=True,
                 # extra arguments
                 **kwargs):
        self.nActors = nActors
        self.percentage_infection = percentage_infection
        self.househouldSize = househouldSize
        self.household_weight = household_weight
        self.household_weight_outbreak = household_weight_outbreak
        self.geo_weight = geo_weight
        self.random_connect_num = random_connect_num
        self.random_weight = random_weight
        self.age_distribution_prob = age_distribution_prob
        self.family_size = family_size
        self.family_age_distribution = family_age_distribution
        self.base_contact_matrix = base_contact_matrix
        self.outbreak_contact_matrix = outbreak_contact_matrix
        self.min_remove = min_remove
        self.household_quar = household_quar
        self.import_cases = import_cases
        self.age_suscep = age_suscep
        self.remove_period = remove_period
        self.outbreak_period = outbreak_period
        self.window_exposed_ratio = window_exposed_ratio
        self.window_incubation = window_incubation
        self.window_post_infectious = window_post_infectious

        self.__dict__.update(kwargs)

    def setparams(self, **kwargs):
        self.__dict__.update(kwargs)

    def simulate(self, seed):
        np.random.seed(2000 + seed)
        acq_network, id_information = acquaintance_network(age_distribution=self.age_distribution_prob,
                                                           base_contact_matrix=self.base_contact_matrix,
                                                           family_size=self.family_size,
                                                           family_age_distribution=self.family_age_distribution,
                                                           nActors=self.nActors,
                                                           householdSize=self.househouldSize,
                                                           weightInhousehold=self.household_weight,
                                                           geo_connection=True,
                                                           geo_weight=self.geo_weight,
                                                           random_connection_num=self.random_connect_num,
                                                           random_weight=self.random_weight)
        day_information, node_information, all_contact_infect \
            = transmission_simulate(acq_network=acq_network,
                                    id_information=id_information,
                                    import_cases=self.import_cases,
                                    base_contact_matrix=self.base_contact_matrix,
                                    outbreak_contact_matrix=self.outbreak_contact_matrix,
                                    age_suscep=self.age_suscep,
                                    remove_period=self.remove_period,
                                    percentage_infection=self.percentage_infection,
                                    weightInhousehold1=self.household_weight,
                                    weightInhousehold2=self.household_weight_outbreak,
                                    household_quar=self.household_quar,
                                    min_remove=self.min_remove,
                                    min_remove2=self.min_remove,
                                    max_day=100,
                                    verbose=False)
        return day_information, node_information, all_contact_infect


class SimulationScenario:

    def __init__(self, scenarios, base_params, n_jobs=16):
        self.scenarios = scenarios
        self.model = AgentBasedSEIR(**base_params)
        self.n_jobs = n_jobs

    def setparams(self, **kwargs):
        self.model.__dict__.update(kwargs)

    def fit(self, i, j):
        start_time = time.time()
        self.model.setparams(**self.scenarios[j])
        result = self.model.simulate(seed=i)
        day_information, id_information, all_contact_infect = result
        day_information['seed'] = i + 1
        id_information['seed'] = i + 1
        all_contact_infect['seed'] = i + 1
        day_information['scenario'] = j
        id_information['scenario'] = j
        all_contact_infect['scenario'] = j
        end_time = time.time()
        results = {"day_information": day_information,
                   "id_information": id_information,
                   "all_contact_infect": all_contact_infect}
        print(f'Scenario {j} Mission {i} completes with running time {(end_time - start_time) / 60: .2f} minutes')
        return results

    def fit_repeat(self, j, n_times):
        with parallel_backend("loky", inner_max_num_threads=1):
            result = Parallel(n_jobs=min(self.n_jobs, n_times), verbose=10)(
                delayed(self.fit)(i, j)
                for i in range(n_times))
        simulate_result = {
            "day_information": pd.concat([x["day_information"] for x in result]),
            "id_information": pd.concat([x["id_information"] for x in result]),
            "all_contact_infect": pd.concat([x["all_contact_infect"] for x in result])
        }
        return simulate_result
