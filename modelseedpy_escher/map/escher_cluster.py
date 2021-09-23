import networkx as nx
import math


class EscherCluster:
    
    def __init__(self, max_distance=50):
        self.max_distance = max_distance
    
    @staticmethod
    def distance(P, Q):
        return math.sqrt(math.pow(Q[0]-P[0], 2) + math.pow(Q[1]-P[1], 2))

    @staticmethod
    def get_cluster(escher_graph, coords_1, max_distance, type_match):
        cluster = set()
        for node_id in escher_graph['nodes']:
            node = escher_graph['nodes'][node_id]
            if node['node_type'] == type_match:
                other_node = node
                coords_2 = (other_node['x'], other_node['y'])
                dist = EscherCluster.distance(coords_1, coords_2)
                if dist < max_distance:
                    cluster.add(node_id)
                    #print(node_id, dist, other_node['bigg_id'])
        return cluster
    
    def compute_all_metabolite_clusters(self, escher_graph, node_type = 'metabolite'):
        g = nx.Graph()
        for node_id in escher_graph['nodes']:
            node = escher_graph['nodes'][node_id]
            if not node_id in g.nodes and node['node_type'] == node_type:
                node = escher_graph['nodes'][node_id]
                coords_1 = (node['x'], node['y'])
                cluster = self.get_cluster(escher_graph, coords_1, self.max_distance, node_type)
                cluster = list(cluster)
                if len(cluster) > 1:
                    prev = cluster[0]
                    for i in range(len(cluster) - 1):
                        g.add_edge(prev, cluster[i + 1])

        cc = list(nx.algorithms.connected_components(g))
        return cc
    
    def cluster(self, escher_map):
        return self.compute_all_metabolite_clusters(escher_map.escher_graph)
    
    def cluster_ids(self, cluster_data, em):
        g = nx.Graph()
        for c in cluster_data:
            if len(c) > 0:
                p = em.nodes[list(c)[0]]['bigg_id']
                for uid in c:
                    id = em.nodes[uid]['bigg_id']
                    #print(p, id)
                    g.add_edge(p, id)
        cc = list(nx.algorithms.connected_components(g))
        return cc
    
    def ids_to_uid(self, cc, em):
        ret = []
        id_to_uid = {}
        for uid in em.nodes:
            n = em.nodes[uid]
            if n['node_type'] == 'metabolite':
                id_to_uid[n['bigg_id']] = uid

        for o in cc:
            c = set()
            for id in o:
                c.add(id_to_uid[id])
            if len(c) > 0:
                ret.append(c)
        return ret
