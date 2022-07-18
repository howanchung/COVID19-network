import itertools
import time
import math
import numpy as np
import numba as nb
import pandas as pd
from scipy.stats import weibull_min


def transmission_simulate(acq_network,
                          id_information,
                          import_cases,
                          base_contact_matrix,
                          outbreak_contact_matrix,
                          age_suscep,
                          remove_period,
                          window_exposed_ratio=None,
                          window_incubation=None,
                          window_post_infectious=None,
                          percentage_infection=0.05,
                          weightInhousehold1=5,
                          weightInhousehold2=10,
                          household_quar=False,
                          min_remove=1,
                          min_remove2=1,
                          outbreak_period=None,
                          max_day=math.inf,
                          pure_infect=True,
                          verbose=True):
    if window_post_infectious is None:
        window_post_infectious = [5.0, 1.4]
        # window_post_infectious = [4, 3.5 / 4]
    if window_incubation is None:
        # code of Temporal dynamics in viral shedding and transmissibility of COVID-19
        # on log-scale
        window_incubation = np.log(np.exp([1.434065, 0.6612]))
    if window_exposed_ratio is None:
        # latent period has mean 4.19 and variance 1.53 according to Zhu et al. (2021)
        # if also assuming a log-normal distribution
        # the ratio between latent and incubation period is also log-normal

        # on log-scale
        window_latent = np.asarray([np.log(4.19) - 0.5 * np.log((1.53 ** 0.5 / 4.19) ** 2 + 1),
                                    np.sqrt(np.log((1.53 ** 0.5 / 4.19) ** 2 + 1))])
        # on log-scale, latent_period / incubation_period
        window_exposed_ratio = [window_incubation[0] - window_latent[0],
                                (window_incubation[1] ** 2 - window_latent[1] ** 2) ** 0.5]
    if outbreak_period is None:
        outbreak_period = [16, 32]

    n = acq_network.shape[0]
    node = np.arange(n)
    # incubation period
    node_incubation = np.minimum(np.random.lognormal(mean=window_incubation[0],
                                                     sigma=window_incubation[1], size=n), 27)
    # latent period
    node_window_exposed = (node_incubation
                           / np.maximum(np.random.lognormal(mean=window_exposed_ratio[0],
                                                            sigma=window_exposed_ratio[1],
                                                            size=n), 1))
    exceed_id = np.where((node_incubation - node_window_exposed) >= 12.3)[0]
    node_window_exposed[exceed_id] = node_incubation[exceed_id] - 12.3
    # pre-symptomatic infectious period
    node_pre_infectious = node_incubation - node_window_exposed
    # post-symptomatic infectious period
    node_infectious_period = np.random.gamma(shape=window_post_infectious[0],
                                             scale=window_post_infectious[1], size=n)
    node_removal_period = np.zeros(n)

    base_contact_matrix = np.maximum(base_contact_matrix, 0.01)
    # node_age = id_information.groupby(by='age')['id'].apply(list)

    infected_day, infectious_day, symp_day, removed_day \
        = np.zeros(n, dtype=int), np.zeros(n, dtype=int), np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    infection_day, iso_day = -np.ones(n), -np.ones(n)

    percentage_infection_node = np.zeros(n)

    all_contact = pd.DataFrame()
    # all_infectious_case = []
    new_symp, count_infections, count_symp = [], [], []
    each_neighbors = acq_network.nonzero()
    each_neighbors = pd.DataFrame({'id': each_neighbors[0], 'contact': each_neighbors[1]})
    each_neighbors = each_neighbors.groupby(by='id')['contact'].apply(np.array)  # equivalent to split in R

    # import cases
    starting_nodes = np.random.choice(node, import_cases.shape[0])
    init_node_index, new_init_node_index = [], []
    import_symp_date = import_cases['onset_date'].to_numpy()
    import_date = np.round(import_symp_date - node_pre_infectious[starting_nodes]).astype(int)
    for i, exposed_day in enumerate(import_date):
        if exposed_day < import_cases.loc[i, 'import_date'] < import_symp_date[i]:
            import_date[i] = import_cases.loc[i, 'import_date']
    day = min(import_date)
    max_import_day = max(import_date)
    all_infected = False

    acq_prob1 = calculate_probas(acq_network.copy())
    acq_prob2 = acq_network.copy()
    acq_prob2[acq_prob2 == weightInhousehold1] = weightInhousehold2
    acq_prob2 = calculate_probas(acq_prob2)

    half_contact_matrix = day_contact(base_contact_matrix, outbreak_contact_matrix, day=10.5)
    min_contact_matrix = outbreak_contact_matrix

    while (not all_infected) and day <= max_day:
        start = time.time()
        if day <= max_import_day:
            new_init_node_index = starting_nodes[np.where(import_date == day)]
            not_yet_index = starting_nodes[np.where(import_date > day)]
            all_index = np.setdiff1d(node, not_yet_index)
            if new_init_node_index.shape[0]:
                infectious_day[new_init_node_index] = np.repeat(day, new_init_node_index.shape[0])
                infected_day[new_init_node_index] = (infectious_day[new_init_node_index] -
                                                     node_window_exposed[new_init_node_index])
                # when imported, cases were already infectious
                infection_day[new_init_node_index] = node_window_exposed[new_init_node_index]
                # infection_day[new_init_node_index] = node_incubation[new_init_node_index]
                node_removal_period[new_init_node_index] = daily_remove(new_init_node_index,
                                                                        day, remove_period, min_remove, min_remove2)
        else:
            all_index = node

        acq_prob = acq_prob1 if day <= outbreak_period[0] or day > outbreak_period[1] else acq_prob2

        if day <= outbreak_period[0]:
            day_contact_matrix = day_contact(base_contact_matrix, outbreak_contact_matrix, day)
        elif day <= outbreak_period[1]:
            day_contact_matrix = min_contact_matrix
        elif day <= outbreak_period[1] + 30:
            day_contact_matrix = day_contact2(half_contact_matrix, min_contact_matrix, 64 - day)
        else:
            day_contact_matrix = half_contact_matrix

        day_contact_num = np.sum(day_contact_matrix, axis=1)
        day_contact_prob = day_contact_matrix / day_contact_num[:, None]  # rowise divide

        all_infectious = np.where(np.logical_and(infection_day >= node_window_exposed,
                                                 infection_day <= node_incubation + node_infectious_period))[0]
        update_infectiousness(percentage_infection_node,
                              percentage_infection, infection_day, node_incubation)
        all_new_infected = np.array([])
        if pure_infect:
            target_member = all_index
        else:
            target_member = node.copy()
        if all_infectious.shape[0]:
            all_contact_day = [individual_contact(node_id=i,
                                                  all_index=target_member,
                                                  node_id_neighbors=each_neighbors[i],
                                                  id_information_age=id_information['age'],
                                                  acq_weight=acq_prob,
                                                  base_contact_num=day_contact_num,
                                                  base_contact_prob=day_contact_prob)
                               for i in all_infectious]
            all_contact_day = pd.concat(all_contact_day, axis=0, ignore_index=True)
            if all_contact_day.shape[0]:
                all_contact_day['suscep'] = age_suscep.loc[
                    all_contact_day['id_age'].to_numpy() - 1, 'suscep'].to_numpy()
                percentage_infection_list = percentage_infection_node[all_contact_day['id'].to_numpy()]
                all_contact_day['infect_prob'] = all_contact_day['suscep'] * percentage_infection_list
                all_contact_day['infect'] = np.concatenate([np.random.choice([True, False], size=1, p=[i, 1 - i])
                                                            for i in all_contact_day['infect_prob']])
                all_contact_day['susceptible'] = (infection_day[all_contact_day['contact']] == -1)
                all_contact_day['contact_freq'] = (all_contact_day
                                                   .groupby(by=['id', 'contact'])['id_age']
                                                   .transform('count'))
                all_contact_day['infect_prob'] = (1 -
                                                  (1 - all_contact_day['infect_prob'])
                                                  ** all_contact_day['contact_freq'].to_numpy())
                all_contact_day = all_contact_day.groupby(by=['id', 'contact']). \
                    agg({'id_age': np.max,
                         'contact_age': np.max,
                         'suscep': np.max,
                         'infect_prob': np.max,
                         'susceptible': np.max,
                         'contact_freq': np.max,
                         'infect': np.any}).reset_index()
                all_contact_day['infect'] = np.logical_and(all_contact_day['infect'],
                                                           all_contact_day['susceptible'])
                all_contact_day['day'] = day
                all_contact = pd.concat([all_contact, all_contact_day], axis=0, ignore_index=True)
                all_new_infected = all_contact_day.loc[all_contact_day['infect'], 'contact'].to_numpy()

        old_infectious = np.where(infection_day >= node_window_exposed)[0]
        old_symp = np.where(infection_day >= node_incubation)[0]
        old_removed = np.where(infection_day >= node_incubation + node_removal_period)[0]
        infection_day[infection_day >= 0] += 1
        if all_new_infected.shape[0]:
            infection_day[all_new_infected] = 0
            infected_day[all_new_infected] = day

        all_infectious_case = np.where(infection_day >= node_window_exposed)[0]
        count_infections.append(all_infectious_case.shape[0])

        all_symp_case = np.where(infection_day >= node_incubation)[0]
        count_symp.append(all_symp_case.shape[0])

        new_infectious_day = list_notcontain(all_infectious_case, old_infectious) \
            if old_infectious.shape[0] else all_infectious_case
        if new_infectious_day.shape[0]:
            infectious_day[new_infectious_day] = day

        new_symp_day = list_notcontain(all_symp_case, old_symp) if old_symp.shape[0] else all_symp_case
        if new_symp_day.shape[0]:
            symp_day[new_symp_day] = day
            node_removal_period[new_symp_day] = daily_remove(new_symp_day, day, remove_period, min_remove,
                                                             min_remove2)
        new_symp.append(np.hstack((new_symp_day, new_init_node_index)))
        new_removed = np.where(infection_day >= node_incubation + node_removal_period)[0]
        if old_removed.shape[0]:
            new_removed = list_notcontain(new_removed, old_removed)

        if household_quar:
            iso_day[iso_day >= 0] += 1
            if new_removed.shape[0]:
                # new_isolated: list = list(set(concatenate_with_none(
                #     [find_family_member(i, id_information) for i in new_removed])))
                new_isolated = np.unique(np.hstack([find_family_member(i, id_information) for i in new_removed]))
                new_isolated_removed = np.asarray([x for x in new_isolated if infection_day[x] > 0])
                if new_isolated_removed.shape[0]:
                    # new_removed = np.array([*new_removed, *new_isolated_removed])
                    new_removed = np.hstack((new_removed, new_isolated_removed))

                if new_isolated.shape[0] and new_isolated_removed.shape[0]:
                    new_isolated = list_notcontain(new_isolated, new_isolated_removed)
                if new_isolated.shape[0]:
                    iso_day[new_isolated] = 0
                    infection_day[new_isolated] = -2
            iso_release = np.where(iso_day == 14)[0]
            if iso_release.shape[0]:
                infection_day[iso_release], iso_day[iso_release] = -1, -1

        if new_removed.shape[0]:
            infection_day[new_removed] = -2
            removed_day[new_removed] = day

        all_infected = all(infection_day < 0)
        end = time.time()
        if verbose:
            print(f'Day= {day}, Active Symp= {count_symp[day - 1]}, Running time= {(end - start): .3f} sec')
        day += 1

    new_symp = [len(x) for x in new_symp]
    day_information = pd.DataFrame({'new_symp': new_symp,
                                    'count_infections': count_infections,
                                    'count_symp': count_symp,
                                    'day': np.arange(len(count_symp)) + min(import_date)})
    node_information = pd.DataFrame({'id': node,
                                     'household': id_information['household'],
                                     'infectious_day': infectious_day,
                                     'symp_day': symp_day,
                                     'infected_day': infected_day,
                                     'removed_day': removed_day})
    sub_columns = ['id', 'contact', 'id_age', 'contact_age', 'contact_freq', 'day']
    all_contact_infect = all_contact.query("infect == True").loc[:, sub_columns].reset_index(drop=True)
    return day_information, node_information, all_contact_infect


def individual_contact(node_id, all_index, node_id_neighbors,
                       id_information_age, acq_weight, base_contact_num, base_contact_prob):
    # neighbors = np.asarray(list_contain(node_id_neighbors, all_index)) if all_index.shape[0] else node_id_neighbors
    neighbors = np.intersect1d(node_id_neighbors, all_index, assume_unique=True) \
        if all_index.shape[0] else node_id_neighbors
    id_acq_weight = acq_weight[node_id, neighbors].toarray().flatten()  # sparse matrix to numpy array
    id_probas = normalize(id_acq_weight)

    id_age = id_information_age[node_id]
    id_neighbors_age = pd.DataFrame({'id_neighbors': neighbors,
                                     'neighbors_age':
                                         id_information_age[neighbors]}).reset_index(drop=True)

    # id_neighbors_age_dist = collections.Counter(id_neighbors_age['neighbors_age'])
    # age_group = np.array(range(1, 15))
    # no_connection_age = list_notcontain(age_group, np.array([*id_neighbors_age_dist.keys()]))
    id_neighbors_age_dist = np.bincount(id_neighbors_age['neighbors_age'], minlength=15)
    no_connection_age = np.where(id_neighbors_age_dist == 0)[0][1:]

    id_contact_num = np.random.poisson(base_contact_num[id_age - 1], size=1)[0]
    id_neighbors_age = id_neighbors_age.groupby("neighbors_age")['id_neighbors'].apply(np.array).to_dict()

    if id_contact_num:
        re_sample = True
        goal_age_count = np.zeros(15, dtype=int)
        while re_sample:
            goal_age_count = np.random.choice(range(1, 15), size=id_contact_num,
                                              replace=True, p=base_contact_prob[id_age - 1, :])
            # goal_age_count = collections.Counter(goal_age_count)
            # re_sample = any([goal_age_count[k] > 0 for k in no_connection_age])
            goal_age_count = np.bincount(goal_age_count, minlength=15)
            re_sample = goal_age_count[no_connection_age].any()

        id_contact = np.random.choice(neighbors, size=id_contact_num, replace=True, p=id_probas)
        id_contact_age = pd.DataFrame({'id': node_id,
                                       'contact': id_contact,
                                       'id_age': id_age,
                                       'contact_age': id_information_age[id_contact]}).reset_index(drop=True)
        # id_contact_age_count = collections.Counter(id_contact_age['contact_age'])
        id_contact_age_count = np.bincount(id_contact_age['contact_age'], minlength=15)
        # age_diff_count = [id_contact_age_count[i] - goal_age_count[i] for i in range(1, 15)]
        age_diff_count = id_contact_age_count - goal_age_count

        add_edges, prune_edges = [], []

        for i in range(1, 15):
            diff_count: int = age_diff_count[i]
            if not diff_count:
                continue
            elif diff_count > 0:
                all_edges_i = np.where(id_contact_age['contact_age'] == i)[0]
                if all_edges_i.shape[0] == 1:
                    prune_edges.append(all_edges_i)
                else:
                    prune_edges.append(np.random.choice(all_edges_i, size=diff_count))
            else:
                id_neighbors_age_i = id_neighbors_age[i]
                if id_neighbors_age_i.shape[0] == 1:
                    add_edges.append(np.repeat(id_neighbors_age_i, -diff_count))
                else:
                    add_age_i_prob = acq_weight[node_id, id_neighbors_age_i].toarray().flatten()
                    # add_age_i_prob = add_age_i_prob / sum(add_age_i_prob)
                    add_age_i_prob = normalize(add_age_i_prob)
                    add_edges.append(np.random.choice(id_neighbors_age_i, size=(-diff_count),
                                                      p=add_age_i_prob, replace=True))

        add_edges = np.concatenate(add_edges) if add_edges else np.array([])
        prune_edges = np.concatenate(prune_edges) if prune_edges else np.array([])

        if prune_edges.shape[0]:
            id_contact_age.drop(index=prune_edges, inplace=True)

        id_contact_age = id_contact_age[['id', 'contact', 'id_age', 'contact_age']]

        if add_edges.shape[0]:
            add_edges_df = pd.DataFrame({'id': node_id,
                                         'contact': add_edges,
                                         'id_age': id_age,
                                         'contact_age': list(id_information_age[add_edges])})
            if not id_contact_age.shape[0]:
                id_contact_age = add_edges_df
            else:
                id_contact_age = pd.concat((id_contact_age, add_edges_df), axis=0, ignore_index=True)
        return id_contact_age
    else:
        return pd.DataFrame()


def daily_contact_matrix(day, outbreak_period, base_contact_matrix, outbreak_contact_matrix,
                         half_contact_matrix, min_contact_matrix=None):
    if min_contact_matrix is None:
        min_contact_matrix = outbreak_contact_matrix
    if day <= outbreak_period[0]:
        day_contact_matrix = day_contact(base_contact_matrix, outbreak_contact_matrix, day)
    elif day <= outbreak_period[1]:
        day_contact_matrix = min_contact_matrix
    elif day <= outbreak_period[1] + 30:
        day_contact_matrix = day_contact2(half_contact_matrix, min_contact_matrix, 64 - day)
    else:
        day_contact_matrix = half_contact_matrix
    return day_contact_matrix


def update_infectiousness(infectiousness_array, percentage_infection, infection_day, node_incubation):
    __curve = {-4: 0.86, -3: 1.03, -2: 1.28, -1: 1.3,
               0: 1.34, 1: 1.33, 2: 1.27,
               3: 1.19, 4: 1.11, 5: 1.05}
    infected = np.where(infection_day > 0)[0]
    for i in infected:
        tmp = round(infection_day[i] - node_incubation[i])
        if tmp <= -5:
            infectiousness = 0.78
        else:
            infectiousness = __curve.get(tmp, 1)
        infectiousness_array[i] = infectiousness * percentage_infection
    return infectiousness_array


@nb.jit(nopython=True, nogil=True)
def day_contact_element(day, c1, c2, d, m):
    c1 = max(c1, 0.01)
    alp = [d, m, c2 / c1, c1]
    return alp[3] * ((1 - alp[2]) / (1 + np.exp(2 * np.log(99) / alp[1] * (day - alp[0] - alp[1] / 2))) + alp[2])


@nb.jit(nopython=True, nogil=True)
def day_contact(base_contact_matrix, outbreak_contact_matrix, day):
    base_contact_matrix_array = base_contact_matrix.flatten()
    outbreak_contact_matrix_array = outbreak_contact_matrix.flatten()
    day_contact_list = [day_contact_element(day, base_contact_matrix_array[i],
                                            outbreak_contact_matrix_array[i],
                                            d=3, m=13) for i in range(14 ** 2)]
    day_contact_matrix = np.array(day_contact_list).reshape(14, 14)
    return day_contact_matrix


@nb.jit(nopython=True, nogil=True)
def day_contact2(base_contact_matrix, outbreak_contact_matrix, day):
    base_contact_matrix_array = base_contact_matrix.flatten()
    outbreak_contact_matrix_array = outbreak_contact_matrix.flatten()
    day_contact_list = [day_contact_element(day, base_contact_matrix_array[i],
                                            outbreak_contact_matrix_array[i],
                                            d=1, m=30) for i in range(14 ** 2)]
    day_contact_matrix = np.array(day_contact_list).reshape(14, 14)
    return day_contact_matrix


def find_family_member(remove_id, id_information):
    family_id = id_information.query(f"id == {remove_id}")['household'].tolist()[0]
    family_member = id_information.query(f"household == {family_id}")['id'].to_numpy()
    return np.setdiff1d(family_member, remove_id)


def daily_remove(new_infection, day, remove_period, min_remove, min_remove2):
    if day <= max(remove_period['date']):
        ref_remove: np.array = remove_period.loc[remove_period['date'] == day, 'time_gap'].values
        remove_p: float = max(ref_remove[0], min_remove) if ref_remove.shape[0] else min_remove
        if not remove_p:
            node_window_recovery = weibull_min.rvs(5, scale=1 / math.gamma(1 + 1 / 5),
                                                   size=new_infection.shape[0])
        else:
            node_window_recovery = weibull_min.rvs(5, scale=remove_p / math.gamma(1 + 1 / 5),
                                                   size=new_infection.shape[0])
    else:
        node_window_recovery = weibull_min.rvs(5, scale=min_remove2 / math.gamma(1 + 1 / 5),
                                               size=new_infection.shape[0])
    return node_window_recovery


def concatenate_with_none(x):
    return np.array(itertools.chain.from_iterable([i for i in x if i]))  # concatenate


def calculate_probas(network):
    network[network > 0] = np.exp(network[network > 0])
    return network


# @nb.jit(nopython=True, nogil=True)
# def list_contain(list1, list2):
#     return [i for i in list1 if i in list2]


@nb.jit(nopython=True, nogil=True)
def list_notcontain(list1, list2):
    return np.asarray([i for i in list1 if i not in list2])


@nb.jit(nopython=True, nogil=True)
def normalize(array):
    return array / np.sum(array)
