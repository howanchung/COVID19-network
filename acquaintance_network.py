import itertools
import numpy as np
import numba as nb
import pandas as pd
from scipy.sparse import csr_matrix

INFLATION_FACTOR = 3


@nb.jit(nopython=True, nogil=True)
def l2_dist(x, y, i):
    return np.sqrt(np.square(x[i] - x) + np.square(y[i] - y))


def age_connect(j, id_age, i_dist, node_age_j, group_num_j):  # j from 1 to 14
    j_dist = node_age_j[np.argsort(i_dist[node_age_j])]
    i_connect_j = j_dist[:round(INFLATION_FACTOR * group_num_j)] if j != id_age \
        else j_dist[1:round(INFLATION_FACTOR * group_num_j + 1)]
    return np.random.choice(i_connect_j, size=group_num_j, replace=False)


def geo_connect(i, id_age, geo_x, geo_y, base_contact_prob, base_contact_num, node_age):
    i_dist = l2_dist(geo_x, geo_y, i)
    i_age = id_age[i]
    base_num, base_prob = base_contact_num[i_age - 1], base_contact_prob[i_age - 1]
    a = np.random.choice(range(1, 15), replace=True, size=round(base_num * INFLATION_FACTOR), p=base_prob)
    # group_num = collections.Counter(a)
    group_num = np.bincount(a, minlength=15)
    i_connect = [age_connect(j, i_age, i_dist, node_age[j - 1], group_num[j]) for j in range(1, 15) if group_num[j]]
    i_connect = list(itertools.chain.from_iterable(i_connect))
    return [i] * len(i_connect), i_connect


def replicate_each(vec, n):
    return list(itertools.chain.from_iterable(itertools.repeat(x, n) for x in vec))


@nb.jit(nopython=True, nogil=True)
def list_notequal(list1, num):
    return [i for i in list1 if i is not num]


@nb.jit(nopython=True, nogil=True)
def multi_sample(pool, col, row):
    samples = []
    for _ in range(row):
        samples.append(np.random.choice(pool, size=col, replace=False))
    return samples


def distribute_age(id_information, family_size, family_age_distribution, age_distribution):
    network_family_size = id_information['household'].value_counts().reset_index()  # equivalent to table in R
    network_family_size.columns = ['family_index', 'family_size']
    network_family_size['reference_index'], network_family_size['family_index2'] = 0, 0
    network_family_size['family_index'] = network_family_size['family_index'].astype(int)
    network_family_size.sort_values(by='family_size', axis=0, inplace=True)
    network_family_size.reset_index(inplace=True)

    unique_observed_size = family_size['family_size'].unique()
    for observed_size in unique_observed_size:
        existing_index = (network_family_size['family_size'] == observed_size)
        observed_size_number = family_size.loc[family_size['family_size'] == observed_size, 'family_index']

        if observed_size_number.shape[0] > 1:
            network_family_size.loc[existing_index, 'reference_index'] = np.random.choice(observed_size_number,
                                                                                          size=np.sum(existing_index))
        else:
            network_family_size.loc[existing_index, 'reference_index'] = [observed_size_number] * np.sum(existing_index)

    network_family_size['reference_index'] = network_family_size['reference_index'].astype(int)

    # not_unique_size = np.setdiff1d(
    #     network_family_size.loc[network_family_size['reference_index'] == 0, 'family_size'].unique(), 1, True)

    not_unique_size = list_notequal(
        network_family_size.loc[network_family_size['reference_index'] == 0, 'family_size'].unique(), 1)

    individual_index = network_family_size.loc[network_family_size['family_size'] == 1, 'family_index']

    for nus in not_unique_size:
        b = nus - unique_observed_size
        b = abs(max(b[b < 0]))
        reference_size = nus + b

        observed_family_index_i = family_size.loc[family_size['family_size'] == reference_size, 'family_index']
        network_family_index_i = network_family_size.loc[
            network_family_size['family_size'] == nus, 'family_index']

        ref_individual = np.random.choice(individual_index, size=b * network_family_index_i.shape[0],
                                          replace=False)
        ref_individual = pd.DataFrame({'ref_id': ref_individual,
                                       'family_id': replicate_each(network_family_index_i, b)})

        reference_i = np.random.choice(observed_family_index_i, size=network_family_index_i.shape[0]) \
            if observed_family_index_i.shape[0] > 1 \
            else [observed_family_index_i] * network_family_index_i.shape[0]

        # reference_i has equal length with network_family_index_i
        network_family_size.loc[network_family_size['family_index'].isin(network_family_index_i),
                                'reference_index'] = reference_i  # set reference index for existing family
        ref_individual['reference_i'] = replicate_each(reference_i, b)  # set reference index for individual
        network_family_size.loc[network_family_size['family_index'].isin(ref_individual['ref_id']),
                                'reference_index'] = ref_individual['reference_i'].to_numpy(dtype='int64')
        network_family_size.loc[network_family_size['family_index'].isin(ref_individual['ref_id']),
                                'family_index2'] = ref_individual['family_id'].to_numpy(dtype='int64')
        individual_index = network_family_size.loc[(network_family_size['family_size'] == 1)
                                                   & (network_family_size['family_index2'] == 0), 'family_index']

    matched_id = (network_family_size['family_index2'] == 0)
    network_family_size.loc[matched_id, 'family_index2'] = network_family_size.loc[matched_id, 'family_index']

    id_information2 = id_information.copy()
    id_information2['household2'] = id_information2['household'].replace(
        network_family_size['family_index'].tolist(),
        network_family_size['family_index2'].tolist())
    # equivalent to mapvalues in R. from value to value
    id_information2['age'] = 0
    unique_household_id = id_information2['household2'].unique()

    for uhi in unique_household_id:
        ref_i = np.unique(network_family_size.loc[network_family_size['family_index2'] == uhi,
                                                  'reference_index'])[0]
        if ref_i:
            a = id_information2.loc[id_information2['household2'] == uhi, 'age'].shape[0]
            b = family_age_distribution.loc[
                family_age_distribution['family_index2'] == ref_i, 'age_group'].shape[0]
            if a != b:
                print(a, b)  # TODO: there is a bug here, the following two array may have different size
                id_information2.loc[id_information2['household2'] == uhi, 'age'] \
                    = family_age_distribution.loc[
                          family_age_distribution[
                              'family_index2'] == ref_i, 'age_group'].to_numpy(dtype='int64')[:a]
            else:
                id_information2.loc[id_information2['household2'] == uhi, 'age'] = family_age_distribution.loc[
                    family_age_distribution['family_index2'] == ref_i, 'age_group'].to_numpy(dtype='int64')

    id_information = id_information2[['id', 'household', 'age']].copy()
    ind_id = (id_information['age'] == 0)
    id_information.loc[ind_id, 'age'] = np.random.choice(range(1, 15),
                                                         size=sum(ind_id),
                                                         replace=True,
                                                         p=age_distribution['freq'].to_numpy())
    return id_information


def acquaintance_network(age_distribution,
                         base_contact_matrix,
                         family_size,
                         family_age_distribution,
                         nActors=10000,
                         householdSize=3,
                         weightInhousehold=5,
                         geo_connection=True,
                         geo_weight=0.5,
                         random_connection_num=1,
                         random_weight=0.3):
    # Use a breakpoint in the code line below to debug your script.
    act_id = range(nActors)
    base_contact_matrix = np.maximum(base_contact_matrix, 0.01)
    base_contact_num = base_contact_matrix.sum(axis=1)
    base_contact_prob = base_contact_matrix / base_contact_num[:, None]

    nhousehold = round(nActors / householdSize)
    household_num = np.random.choice(range(nhousehold), size=nActors)

    id_information = pd.DataFrame({'id': act_id, 'household': household_num})

    household_weight = id_information.merge(id_information, left_on='household',
                                            right_on='household',
                                            how='outer')
    household_weight = csr_matrix((np.array([weightInhousehold] * household_weight.shape[0]),
                                   (household_weight['id_x'], household_weight['id_y'])),
                                  shape=(nActors, nActors), dtype=int)
    household_weight.setdiag(0)
    acq_network = household_weight

    family_size = family_size.sort_values(by='family_size', axis=0)
    family_size = family_size.reset_index(drop=True)
    id_information = distribute_age(id_information,
                                    family_size, family_age_distribution, age_distribution)
    node_age = np.array(id_information.groupby(by='age')['id'].apply(np.array))  # equivalent to split in R

    if geo_connection:
        id_information['geo_x'], id_information['geo_y'] = [np.random.uniform(low=0, high=100, size=nActors),
                                                            np.random.uniform(low=0, high=100, size=nActors)]
        geo_connection = [geo_connect(i, id_information['age'].to_numpy(),
                                      id_information['geo_x'].to_numpy(),
                                      id_information['geo_y'].to_numpy(),
                                      base_contact_prob,
                                      base_contact_num, node_age) for i in act_id]
        contact = itertools.chain.from_iterable([i[0] for i in geo_connection])  # concat list of lists
        contacted = itertools.chain.from_iterable([i[1] for i in geo_connection])  # concat list of lists
        geo_connection = pd.DataFrame({'id': contact,
                                       'contact': contacted})
        geo_connection_weight = csr_matrix((np.array([geo_weight] * geo_connection.shape[0]),
                                            (geo_connection['id'], geo_connection['contact'])),
                                           shape=(nActors, nActors), dtype=float)
        geo_connection_weight = geo_connection_weight + geo_connection_weight.T
        acq_network += geo_connection_weight

    if random_connection_num:
        # random_connection = np.array([np.random.choice(range(nActors), size=2, replace=False)
        #                               for _ in range(round(nActors * random_connection_num))])
        random_connection = np.array(multi_sample(np.array(range(nActors)), 2, round(nActors * random_connection_num)))
        random_matrix = csr_matrix((np.array([random_weight] * random_connection.shape[0]),
                                    (random_connection[:, 0], random_connection[:, 1])),
                                   shape=(nActors, nActors), dtype=float)
        random_matrix = random_matrix + random_matrix.T

        acq_network += random_matrix

    acq_network[acq_network >= weightInhousehold] = weightInhousehold

    return acq_network, id_information
