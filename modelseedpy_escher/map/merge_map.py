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
        self.map_ppp = {}
        self.uid_mapping = {}
        self.nodes = {}

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
            logger.debug("%s", em.escher_graph['canvas'])
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

    def compute_segments(self, map_index, reaction, visited, placed):
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

        logger.debug("RXN  %s", reaction['bigg_id'])
        segments = {}
        for seg_uid in reaction['segments']:
            seg = reaction['segments'][seg_uid]
            index_from = self.map_ppp[(map_index, seg['from_node_id'])]
            index_to = self.map_ppp[(map_index, seg['to_node_id'])]
            t = tuple(sorted([min(f(index_from)), min(f(index_to))]))
            logger.debug("SEG  MAP[%d] %s[%d] -> %s[%d] t:%s", map_index,
                           seg['from_node_id'], index_from, seg['to_node_id'], index_to, t)
            if t not in visited:
                logger.debug("ADD  SEG %s %s", index_from, index_to)
                visited.add(t)
                if index_from not in self.uid_mapping:
                    index_from = list(f(index_from) & set(self.uid_mapping))[0]
                if index_to not in self.uid_mapping:
                    index_to = list(f(index_to) & set(self.uid_mapping))[0]
                placed |= {index_from, index_to}
                uid_from = self.uid_mapping[index_from]
                uid_to = self.uid_mapping[index_to]
                segments[str(self.em_uid)] = {
                    'from_node_id': str(uid_from),
                    'to_node_id': str(uid_to),
                    'b1': seg['b1'],
                    'b2': seg['b2']
                }
                self.em_uid += 1
            else:
                logger.debug("SKIP SEG %s %s", index_from, index_to)
        return segments

    def compute_reactions(self, em, map_index, visited, placed):
        reactions = {}
        for reaction in em.reactions:
            # print(reaction['bigg_id'])
            r = copy.deepcopy(reaction)
            segments = self.compute_segments(map_index, reaction, visited, placed)
            if len(segments) > 0:
                r['segments'] = segments
                reactions[str(self.em_uid)] = r
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

    def get_map_stack(self, maps):
        index = 0
        map_stack = []
        for map_index in range(len(maps)):
            em = maps[map_index]
            indexes = set()
            for n in em.nodes:
                indexes.add(index)
                self.nodes[index] = n
                self.map_ppp[(map_index, n['uid'])] = index
                index += 1
            map_stack.append(indexes)
        return map_stack

    def merge(self, maps, max_distance=10):
        self.em_uid = 0
        self.map_ppp = {}
        self.nodes = {}
        self.uid_mapping = {}

        map_stack = self.get_map_stack(maps)

        cc = self.compute_clusters(self.nodes, max_distance)
        self.pointers = self.get_pointers(cc)

        em_result = EscherMap([
            {'map_name': '', 'map_id': '', 'map_description': '', 'homepage': '', 'schema': ''},
            {'reactions': {}, 'nodes': {}, 'text_labels': {}, 'canvas': self.compute_canvas(maps)}
        ])
        visited = set()

        for stack in map_stack:
            for index in stack:
                if index not in visited:
                    # print('copy', nodes[index]['node_type'])
                    em_result.escher_graph['nodes'][self.em_uid] = copy.deepcopy(self.nodes[index])
                    em_result.escher_graph['nodes'][self.em_uid]['uid'] = self.em_uid
                    self.uid_mapping[index] = self.em_uid
                    self.em_uid += 1
                    visited.add(index)
                    if index in self.pointers:
                        visited |= self.pointers[index]
            # copy seg
        logger.debug("Nodes: %d", len(em_result.escher_graph['nodes']))

        reactions = {}
        visited = set()
        placed = set()
        reactions.update(self.compute_reactions(maps[0], 0, visited, placed))
        for map_index in range(1, len(maps)):
            reactions.update(self.compute_reactions(maps[map_index], map_index, visited, placed))

        em_result.escher_graph['reactions'] = reactions

        logger.debug("Reactions: %d", len(em_result.escher_graph['reactions']))
        return em_result
