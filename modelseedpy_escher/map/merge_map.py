import logging
import copy
import networkx as nx
from modelseedpy_escher.core import EscherMap
from modelseedpy_escher.map.escher_cluster import EscherCluster

logger = logging.getLogger(__name__)


class EscherMerge:

    def __init__(self):
        self.em_uid = 0
        self.pointers = {}
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

    @staticmethod
    def compute_clusters(nodes, max_distance):
        g = nx.Graph()
        for index in nodes:
            # index_map[index] = set()
            node = nodes[index]
            coords_1 = (node['x'], node['y'])
            cluster = list(EscherMerge.get_cluster(nodes, coords_1, max_distance))
            if len(cluster) > 1:
                prev = cluster[0]
                for i in range(len(cluster) - 1):
                    g.add_edge(prev, cluster[i + 1])
        return list(nx.algorithms.connected_components(g))

    @staticmethod
    def compute_canvas(maps):
        xs = []
        ys = []
        for em in maps:
            print(em.escher_graph['canvas'])
            x = em.escher_graph['canvas']['x']
            y = em.escher_graph['canvas']['y']
            w = em.escher_graph['canvas']['width']
            h = em.escher_graph['canvas']['height']
            xs.append(x)
            xs.append(x + w)
            ys.append(y)
            ys.append(y + h)
        x_min = min(xs)
        x_max = max(xs)
        y_min = min(ys)
        y_max = max(ys)
        return {'x': x_min, 'y': y_min, 'width': x_max - x_min, 'height': y_max - y_min}

    def compute_segments(self, reaction, map_ppp, uid_mapping, visited, placed):
        def f(i):
            """
            Expand to the cluster index set
            return set with i + its pointers

            :param i: index to expand
            :return: set of cluster of indexes
            """
            if i in self.pointers:
                return {i} | self.pointers[i]
            return {i}

        segments = {}
        for seg_uid in reaction['segments']:
            seg = reaction['segments'][seg_uid]
            index_from = map_ppp[(0, seg['from_node_id'])]
            index_to = map_ppp[(0, seg['to_node_id'])]
            t = tuple(sorted([min(f(index_from)), min(f(index_to))]))
            if t not in visited:
                logger.debug("ADD  SEG %s %s", index_from, index_to)
                visited.add(t)
                placed |= {index_from, index_to}
                uid_from = uid_mapping[index_from]
                uid_to = uid_mapping[index_to]
                segments[self.em_uid] = {
                    'from_node_id': uid_from,
                    'to_node_id': uid_to,
                    'b1': seg['b1'],
                    'b2': seg['b2']
                }
                self.em_uid += 1
            else:
                logger.debug("SKIP SEG %s %s", index_from, index_to)
        return segments

    def compute_reactions(self, em, map_ppp, uid_mapping, visited, placed):
        reactions = {}
        for reaction in em.reactions:
            # print(reaction['bigg_id'])
            r = copy.deepcopy(reaction)
            segments = self.compute_segments(reaction, map_ppp, uid_mapping, visited, placed)
            if len(segments) > 0:
                r['segments'] = segments
                reactions[self.em_uid] = r
                self.em_uid += 1
        return reactions

    @staticmethod
    def get_pointers(cc):
        pointers = {}
        for c in cc:
            for index1 in c:
                pointers[index1] = set()
                for index2 in c:
                    if index1 != index2:
                        pointers[index1].add(index2)
        return pointers

    def merge(self, maps, max_distance=10):
        self.em_uid = 0
        nodes = {}
        index = 0
        map_stack = []
        map_ppp = {}
        for map_index in range(len(maps)):
            em = maps[map_index]
            indexes = set()
            for n in em.nodes:
                indexes.add(index)
                nodes[index] = n
                map_ppp[(map_index, n['uid'])] = index
                index += 1
            map_stack.append(indexes)

        cc = self.compute_clusters(nodes, max_distance)
        self.pointers = self.get_pointers(cc)

        em_result = EscherMap([
            {'map_name': '', 'map_id': '', 'map_description': '', 'homepage': '', 'schema': ''},
            {'reactions': {}, 'nodes': {}, 'text_labels': {}, 'canvas': self.compute_canvas(maps)}
        ])
        visited = set()
        import copy
        uid_mapping = {}
        for stack in map_stack:
            for index in stack:
                if index not in visited:
                    # print('copy', nodes[index]['node_type'])
                    em_result.escher_graph['nodes'][self.em_uid] = copy.deepcopy(nodes[index])
                    em_result.escher_graph['nodes'][self.em_uid]['uid'] = self.em_uid
                    uid_mapping[index] = self.em_uid
                    self.em_uid += 1
                    visited.add(index)
                    if index in self.pointers:
                        visited |= self.pointers[index]
            # copy seg
        logger.debug("Nodes: %d", len(em_result.escher_graph['nodes']))

        reactions = {}
        visited = set()
        placed = set()
        reactions.update(self.compute_reactions(maps[0], map_ppp, uid_mapping, visited, placed))
        for map_index in range(1, len(maps)):
            reactions.update(self.compute_reactions(maps[map_index], map_ppp, uid_mapping, visited, placed))

        em_result.escher_graph['reactions'] = reactions

        logger.debug("Reactions: %d", len(em_result.escher_graph['reactions']))
        return em_result
