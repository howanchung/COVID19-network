import time
from collections import deque
import igraph as ig
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from joblib import Parallel, delayed, parallel_backend


class NetworkAnalyser:

    def __init__(self, result, n_actors=20000, day_limit=60):
        self.result = result
        self.n_actors = n_actors
        self.day_limit = day_limit

    def analyse(self, n_times, n_jobs):
        start_t = time.time()
        with parallel_backend("loky", inner_max_num_threads=1):
            simulation_summary = Parallel(n_jobs=n_jobs, verbose=10)(
                delayed(self.summary_graph_period)
                (i, self.result, self.n_actors, True)
                for i in range(1, n_times + 1))
        simulation_summary_graph_total = pd.concat(simulation_summary)
        simulation_summary_graph_total.reset_index(drop=True, inplace=True)
        simulation_summary_graph_total = self.summary_graph_t_dt(simulation_summary_graph_total)
        end_t = time.time()
        print(f'Finished analyzing network with running time {(end_t - start_t) / 60: .4f} minutes')
        return simulation_summary_graph_total

    @staticmethod
    def summary_graph_t_dt(simulation_summary_dt):
        period = np.sort(simulation_summary_dt['period'].unique())
        a = deque()
        for i in period:
            a_i = simulation_summary_dt.loc[
                simulation_summary_dt['period'] == i, simulation_summary_dt.columns != 'period']
            a_i = a_i.apply(np.quantile, axis=0, q=[0.025, 0.5, 0.975]).T
            a_i.reset_index(inplace=True)
            a_i.insert(0, 'day', i)
            a_i.columns = ['day', 'Measures', 'min', 'median', 'max']
            a.append(a_i)
        a = pd.concat(a)
        a = a.reset_index(drop=True)
        return a

    def summary_graph_period(self, seed, result, n, include_singletons):
        start_time = time.time()
        id_information = result["id_information"].query(f"seed == {seed}")
        day_information = result["day_information"].query(f"seed == {seed}")
        all_contact = result["all_contact_infect"].query(f"seed == {seed}")
        sum_dt_list = [self.sum_graph_dt_period(i, 1, all_contact, day_information,
                                                id_information, n, include_singletons)
                       for i in range(1, self.day_limit + 1)]
        sum_dt_list = pd.concat(sum_dt_list)
        sum_dt_list.fillna(value=0, inplace=True)
        end_time = time.time()
        print(f'Seed {seed} completed with running time {(end_time - start_time) / 60: .4f} minuites')
        return sum_dt_list

    def sum_graph_dt_period(self, time_end, time_begin, all_contact_seed,
                            day_information_seed, id_information_seed, n, include_singletons=False):
        all_contact_interval = (all_contact_seed
                                .query(f"{time_begin} <= day <= {time_end}")
                                .drop_duplicates(subset=['id', 'contact'], ignore_index=True))
        if rows := all_contact_interval.shape[0]:
            # household proportion
            all_contact_interval['id_household'] = \
                (all_contact_interval.loc[:, ['id']]
                 .merge(id_information_seed.loc[:, ["id", "household"]],
                        on='id', how='left'))['household']
            all_contact_interval['contact_household'] = \
                (all_contact_interval.loc[:, ['contact']]
                 .merge(id_information_seed.loc[:, ["id", "household"]],
                        left_on='contact', right_on="id", how='left'))['household']

            household_prop = np.mean(all_contact_interval['id_household'] == all_contact_interval['contact_household'])

            # building networks
            adjacency_matrix = csr_matrix((np.ones(rows, dtype=int),
                                           (all_contact_interval['id'], all_contact_interval['contact'])),
                                          shape=(n, n), dtype=int)
            # only allow one source case
            in_degree = adjacency_matrix.sum(axis=0).A[0]
            node_w_more_than_one_parent = np.where(in_degree > 1)[0]
            if node_w_more_than_one_parent.shape[0]:
                adjacency_matrix = adjacency_matrix.tolil()
                for j in node_w_more_than_one_parent:
                    col_j = adjacency_matrix[:, j].A
                    parent_j = np.where(col_j == 1)[0]
                    true_parent_j = np.random.choice(parent_j, size=1)
                    not_parent_j = np.setdiff1d(parent_j, true_parent_j)
                    adjacency_matrix[not_parent_j, j] = 0
                    adjacency_matrix[true_parent_j, j] = 1
                adjacency_matrix = adjacency_matrix.tocsr()
            all_degree = adjacency_matrix.sum(axis=0).A[0] + adjacency_matrix.sum(axis=1).A.flatten()
            # network attributes
            infected = np.where(
                np.logical_and(
                    id_information_seed['infected_day'] >= time_begin,
                    id_information_seed['infected_day'] <= time_end))[0]
            infected_count = infected.shape[0]

            if time_end <= day_information_seed['new_symp'].shape[0]:
                new_symp = day_information_seed['new_symp'].iloc[time_end - 1]
            else:
                new_symp = 0
            results = {"period": time_end, "Percentage of infection": infected_count / n,
                       "Number of the new-onset": new_symp, "Proportion within household": household_prop}
            results.update(self.get_net_attribute(adjacency_matrix, directed=True))

            if include_singletons:
                singletons = np.where(all_degree == 0)[0]
                singletons = infected[np.isin(infected, singletons)]
                num_singletons = singletons.shape[0]
                p_singletons = np.round(num_singletons / infected_count, 5)
                results.update({"Number of singletons": num_singletons,
                                "Percentage of singletons": p_singletons})
            summary_dt = pd.DataFrame(results, index=[0])
            return summary_dt
        else:
            return pd.DataFrame()

    @staticmethod
    def get_net_attribute(adjacency, directed=True, attribute=None):
        if attribute is None:
            attribute = {"mean_outdegree": True,
                         "dist_T": True,
                         "mean_betweenness": True,
                         "mean_diameter": True,
                         "mean_cluster_size": True}

        mode = 'undirected' if not directed else 'directed'
        graph = ig.Graph.Adjacency(adjacency, mode=mode)
        if graph:
            in_degree, out_degree = adjacency.sum(axis=0).A[0], adjacency.sum(axis=1).A.flatten()
            all_degree = in_degree + out_degree
            non_singletons = np.where(all_degree != 0)
            if attribute['mean_outdegree']:
                mean_outdegree = np.mean(out_degree[non_singletons])
            else:
                mean_outdegree = None
            if attribute['dist_T']:
                dist_T = ig.Graph.average_path_length(graph, directed)
            else:
                dist_T = None
            if any((attribute.get(x) for x in ('mean_diameter', 'mean_cluster_size'))):
                subgraphs = graph.decompose(minelements=2, mode='weak')
                mean_betweenness = np.mean(np.concatenate([x.betweenness() for x in subgraphs]))
                mean_diameter = np.mean([x.diameter(directed) for x in subgraphs])
                mean_cluster_size = np.mean([x.vcount() for x in subgraphs])
            else:
                mean_betweenness = np.mean(graph.betweenness(vertices=non_singletons))
                mean_diameter = graph.diameter(directed)
                mean_cluster_size = None
        else:
            mean_outdegree, dist_T, mean_betweenness, mean_diameter, mean_cluster_size = 0, 0, 0, 0, 0
        results = {"Average outdegree": mean_outdegree,
                   "Average shortest path length": dist_T,
                   "Average betweenness": mean_betweenness,
                   "Average diameter of clusters": mean_diameter,
                   "Average size of clusters": mean_cluster_size}
        return results
