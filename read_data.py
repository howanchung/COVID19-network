import pandas as pd
from daily_contact_network import day_contact, daily_contact_matrix


def read_data():
    import_cases = pd.read_csv('data/import_cases.csv')
    import_cases.drop(import_cases.columns[0], axis=1, inplace=True)
    age_suscep = pd.read_csv('data/age_suscep.csv')
    age_suscep.drop(age_suscep.columns[0], axis=1, inplace=True)
    remove_period = pd.read_csv('data/remove_period.csv')
    remove_period.drop(remove_period.columns[0], axis=1, inplace=True)
    remove_period = remove_period.rename_axis("date").reset_index()
    age_distribution = pd.read_csv('data/ZJ_pop.csv')
    age_distribution_prob = age_distribution.copy()
    age_distribution_prob['freq'] = age_distribution_prob['freq'] / sum(age_distribution_prob['freq'])
    shanghai_contact_base = pd.read_csv('data/shanghai_contact_base.csv')
    shanghai_contact_base.drop(shanghai_contact_base.columns[[0]], axis=1, inplace=True)
    shanghai_contact_base = shanghai_contact_base.to_numpy()
    shanghai_contact_outbreak = pd.read_csv('data/shanghai_contact_outbreak.csv')
    shanghai_contact_outbreak.drop(shanghai_contact_outbreak.columns[[0]], axis=1, inplace=True)
    shanghai_contact_outbreak = shanghai_contact_outbreak.to_numpy()
    family_age_distribution = pd.read_csv('data/family_age_distribution.csv')
    family_age_distribution.drop(family_age_distribution.columns[[0]], axis=1, inplace=True)
    family_size = pd.read_csv('data/family_size.csv')
    family_size.drop(family_size.columns[[0]], axis=1, inplace=True)

    medium_contact_matrix = day_contact(shanghai_contact_base, shanghai_contact_outbreak, day=9.5)
    medium2_contact_matrix = day_contact(medium_contact_matrix, shanghai_contact_outbreak, day=9.5)

    scenarios = {1: {"base_contact_matrix": shanghai_contact_base,
                     "outbreak_contact_matrix": shanghai_contact_outbreak,
                     "min_remove": 1,
                     "household_quar": True},
                 2: {"base_contact_matrix": shanghai_contact_base,
                     "outbreak_contact_matrix": medium_contact_matrix,
                     "min_remove": 3,
                     "household_quar": False},
                 3: {"base_contact_matrix": shanghai_contact_base,
                     "outbreak_contact_matrix": medium_contact_matrix,
                     "min_remove": 3,
                     "household_quar": True},
                 4: {"base_contact_matrix": medium_contact_matrix,
                     "outbreak_contact_matrix": medium_contact_matrix,
                     "min_remove": 3,
                     "household_quar": False},
                 5: {"base_contact_matrix": medium_contact_matrix,
                     "outbreak_contact_matrix": medium2_contact_matrix,
                     "min_remove": 3,
                     "household_quar": False},
                 6: {"base_contact_matrix": medium_contact_matrix,
                     "outbreak_contact_matrix": medium_contact_matrix,
                     "min_remove": 3,
                     "household_quar": True},
                 7: {"base_contact_matrix": medium_contact_matrix,
                     "outbreak_contact_matrix": medium2_contact_matrix,
                     "min_remove": 3,
                     "household_quar": True}}

    base_params = {"nActors": 20000,
                   "percentage_infection": 0.05 / 1.34,
                   "househouldSize": 2.94,
                   "household_weight": 3.716,
                   "household_weight_outbreak": 5.041,
                   "geo_weight": 0.5,
                   "random_weight": 0.5,
                   "random_connect_num": 1,
                   "age_distribution_prob": age_distribution_prob,
                   "family_size": family_size,
                   "family_age_distribution": family_age_distribution,
                   "import_cases": import_cases,
                   "age_suscep": age_suscep,
                   "remove_period": remove_period,
                   "outbreak_period": [16, 32]
                   }

    contact_matrix = [daily_contact_matrix(day, outbreak_period=[16, 32],
                                           base_contact_matrix=shanghai_contact_base,
                                           outbreak_contact_matrix=shanghai_contact_outbreak,
                                           half_contact_matrix=medium_contact_matrix) for day in range(1, 70)]
    return scenarios, base_params, contact_matrix
