import networkx as nx
from modelseedpy_escher.core import EscherMap
from modelseedpy_escher.map.escher_cluster import EscherCluster


class EscherMerge:

    def __init__(self):
        pass

    @staticmethod
    def get_cluster(nodes, coords_1, max_distance):
        cluster = set()
        for index in nodes:
            node = nodes[index]
            coords_2 = (node['x'], node['y'])
            dist = EscherCluster.distance(coords_1, coords_2)
            if dist < max_distance:
                cluster.add(index)
        return cluster

    def merge(self, maps, max_distance=10):
        nodes = {}
        index = 0
        map_stack = []
        for em in maps:
            indexes = set()
            for n in em.nodes:
                indexes.add(index)
                nodes[index] = n
                index += 1
            map_stack.append(indexes)

        g = nx.Graph()
        for index in nodes:
            #index_map[index] = set()
            node = nodes[index]
            coords_1 = (node['x'], node['y'])
            cluster = list(self.get_cluster(nodes, coords_1, max_distance))
            if len(cluster) > 1:
                prev = cluster[0]
                for i in range(len(cluster) - 1):
                    g.add_edge(prev, cluster[i + 1])
        cc = list(nx.algorithms.connected_components(g))

        pointers = {}
        for c in cc:
            for index1 in c:
                pointers[index1] = set()
                for index2 in c:
                    if index1 != index2:
                        pointers[index1].add(index2)

        em3 = EscherMap(
            [{'map_name': '', 'map_id': '', 'map_description': '', 'homepage': '', 'schema': ''},
             {'reactions': {}, 'nodes': {}, 'text_labels': {},
              'canvas': maps[0].escher_graph['canvas']}
             ])
        visited = set()
        em_uid = 0
        import copy
        for stack in map_stack:
            for index in stack:
                if index not in visited:
                    # print('copy', nodes[index]['node_type'])
                    em3.escher_graph['nodes'][em_uid] = copy.deepcopy(nodes[index])
                    em3.escher_graph['nodes'][em_uid]['uid'] = em_uid
                    em_uid += 1
                    visited.add(index)
                    if index in pointers:
                        visited |= pointers[index]
            # copy seg
        print('merge nodes', len(em3.escher_graph['nodes']))

        return em3
