# import os
import time
import numpy as np
import pandas as pd
from read_data import read_data
from ag_based_seir import SimulationScenario
from network_analyser import NetworkAnalyser

N_JOBS = 40
N_TIMES = 200
N_SCENARIOS = 7
PATH = ""

if __name__ == '__main__':
    # os.chdir('C:\\Users\\105c\\Desktop\\python project')
    start_t = time.time()

    scenarios, base_params, contact_matrix = read_data()
    simulation = SimulationScenario(scenarios=scenarios, base_params=base_params, n_jobs=N_JOBS)

    simulate_result = {}
    for j in range(1, N_SCENARIOS + 1):
        simulate_result[j] = simulation.fit_repeat(j, n_times=N_TIMES)

    contact_pattern = pd.DataFrame(np.vstack(np.sum(contact_matrix, axis=1)))
    contact_pattern.to_csv(PATH + "contact_pattern.csv", index=False, header=False)

    day_information = pd.concat([x["day_information"] for _, x in simulate_result.items()])
    id_information = pd.concat([x["id_information"] for _, x in simulate_result.items()])
    all_contact_infect = pd.concat([x["all_contact_infect"] for _, x in simulate_result.items()])

    day_information.to_parquet(PATH + 'day_information.parquet')
    id_information.to_parquet(PATH + 'id_information.parquet')
    all_contact_infect.to_parquet(PATH + 'all_contact_infect.parquet')
    end_t = time.time()
    print(f'Total running time: {(end_t - start_t) / 3600: .4f} hours')

    result = {"day_information": pd.read_parquet(PATH + 'day_information.parquet'),
              "id_information": pd.read_parquet(PATH + 'id_information.parquet'),
              "all_contact_infect": pd.read_parquet(PATH + 'all_contact_infect.parquet')}

    simulation_summary_total = []
    for j in range(1, N_SCENARIOS + 1):
        print(f'Starting analyzing scenario {j}')
        start_t = time.time()
        result_j = {k: v.query(f"scenario == {j}") for k, v in result.items()}
        analyser = NetworkAnalyser(result_j, day_limit=100)
        simulation_summary = analyser.analyse(n_times=N_TIMES, n_jobs=N_JOBS)
        simulation_summary['scenario'] = j
        simulation_summary_total.append(simulation_summary)
        end_t = time.time()
        print(f'Finished analyzing scenario {j} with running time {(end_t - start_t) / 3600: .4f} hours')

    simulation_summary_graph_total = pd.concat(simulation_summary_total)
    simulation_summary_graph_total.to_csv(PATH + 'simulation_summary_total.csv', index=False)
