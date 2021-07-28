import logging
import copy

logger = logging.getLogger(__name__)

class EscherGrid:
    
    def __init__(self):
        pass
    
    def get_next_id(self, em):
        next_id = 0
        for node_id in em.escher_graph['nodes']:
            if int(node_id) >= next_id:
                next_id = int(node_id) + 1

        return next_id
    
    def copy_map_compounds(self, escher_map, escher_map_target, next_id, cpd_lookup, x_off = 0, y_off = 0):
        compound_remap = {}
        node_id_remap = {}
        id_counter = next_id
        for node_id in escher_map[1]['nodes']:
            node = escher_map[1]['nodes'][node_id]
            note_type = node['node_type']
            if note_type == 'metabolite':
                #print(node)
                tnode = self.translate_metabolite_node(node, cpd_lookup, x_off, y_off)
                #print(node, tnode)
                #print(node['bigg_id'], tnode['bigg_id'])
                compound_remap[node['bigg_id']] = tnode['bigg_id']
                escher_map_target[1]['nodes'][str(id_counter)] = tnode
            else:
                escher_map_target[1]['nodes'][str(id_counter)] = {
                    'x': node['x'] + x_off,
                    'y': node['y'] + y_off,
                    'node_type': node['node_type']
                }
            node_id_remap[node_id] = str(id_counter)
            id_counter += 1

        return id_counter, node_id_remap, compound_remap
    
    def translate_metabolite_node(self, n, lookup, x_off = 0, y_off = 0):
        #print(n)
        bigg_id, name = lookup(n)
        #print(bigg_id, name)
        #id = node['bigg_id']
        #id = id[:-2]
        if bigg_id == None:
            raise Exception('not found ' + n['bigg_id'])

        result = copy.deepcopy(n)
        result['x'] += x_off
        result['y'] += y_off
        result['label_x'] += x_off
        result['label_y'] += y_off
        result['bigg_id'] = bigg_id
        result['name'] = name
        #result = {
        #    'x': n['x'] + x_off, 
        #    'y': n['y'] + y_off, 
        #    'name': name, 
        #    'node_is_primary': n['node_is_primary'], 
        #    'node_type': 'metabolite', 
        #    'label_x': n['label_x'] + x_off, 
        #    'label_y': n['label_y'] + y_off,
        #    'bigg_id': bigg_id
        #}
        return result

    def translate_reaction_node(self, rnode, lookup, next_id, x_off = 0, y_off = 0,
                                node_id_remap = {}, compound_remap = {}, allow_miss = True):
        #print(compound_remap)
        id = rnode['bigg_id']
        bigg_id, name = lookup(rnode)
        if bigg_id == None:
            if allow_miss:
                bigg_id = rnode['bigg_id']
                name = rnode['name']
            else:
                raise Exception('not found ' + rnode['bigg_id'])

        metabolites = []
        segments = {}
        for m in rnode['metabolites']:
            cpd_id = compound_remap[m['bigg_id']] if m['bigg_id'] in compound_remap else '?'
            metabolites.append({
                'bigg_id' : cpd_id,
                'coefficient' : m['coefficient']
            })
        for s in rnode['segments']:
            segment = rnode['segments'][s]
            if segment['from_node_id'] in node_id_remap and segment['to_node_id'] in node_id_remap:
                segment_translate = {
                    'b1' : self.translate_b(segment['b1'], x_off, y_off),
                    'b2' : self.translate_b(segment['b2'], x_off, y_off),
                    'from_node_id' : node_id_remap[segment['from_node_id']],
                    'to_node_id' : node_id_remap[segment['to_node_id']]
                }
                segments[str(next_id)] = segment_translate
                next_id += 1
        result = {
            'bigg_id': bigg_id,
            'name': name,
            'gene_reaction_rule': '',
            'genes': [],
            'label_x': rnode['label_x'] + x_off,
            'label_y': rnode['label_y'] + y_off,
            'metabolites': metabolites,
            'reversibility': rnode['reversibility'],
            'segments': segments
        }
        return result, next_id
    
    def copy_map_reactions(self, escher_map, escher_map_target, next_id, node_id_remap, compound_remap, rxn_lookup, x_off = 0, y_off = 0):
        id_counter = next_id
        for rnode_id in escher_map[1]['reactions']:
            rnode = escher_map[1]['reactions'][rnode_id]

            trnode, id_counter = self.translate_reaction_node(rnode, rxn_lookup, id_counter, 
                                                      x_off, y_off, 
                                                      node_id_remap, compound_remap)

            escher_map_target[1]['reactions'][str(id_counter)] = trnode
            id_counter += 1

        return id_counter

    def translate_b(self, segment, x_off, y_off):
        if segment == None:
            return None
        else:
            return {
                'x' : segment['x'] + x_off, 
                'y' : segment['y'] + y_off}
    
    @staticmethod
    def get_max_canvas(maps):
        max_w = 0
        max_h = 0
        for e in maps:
            w = e.escher_graph['canvas']['width']
            h = e.escher_graph['canvas']['height']
            if w > max_w:
                max_w = w
            if h > max_h:
                max_h = h
        return max_w, max_h
    
    def build(self, maps, grid):
        master = maps[0].clone()
        master.set_to_origin()
        block_width, block_height = self.get_max_canvas(maps)
        logger.debug('Max W: %f, H: %f', block_width, block_height)
        master.escher_graph['canvas']['width'] = block_width * grid[0]
        master.escher_graph['canvas']['height'] = block_height * grid[1]
        master.escher_graph['canvas']['x'] = 0
        master.escher_graph['canvas']['y'] = 0
        #block_width = master.escher_graph['canvas']['width']
        #block_height = master.escher_graph['canvas']['height']
        #master.escher_graph['canvas']['width'] *= grid[0]
        #master.escher_graph['canvas']['height'] *= grid[1]
        
        next_id = self.get_next_id(maps[0])
        
        logger.debug('next_id: %d', next_id)

        map_index = 1
        for i in range(0, grid[1]):
            for j in range(0, grid[0]):
                logger.debug("[%d, %d]", i, j)
                if not (i == 0 and j == 0) and map_index < len(maps):
                    if not maps[map_index] == None:
                        to_draw = maps[map_index].clone()
                        to_draw.set_to_origin()
                        logger.debug("[DRAW] %d next_id: %d", map_index, next_id)

                        next_id, node_id_remap, compound_remap = self.copy_map_compounds(to_draw.escher_map,
                                           master.escher_map, 
                                           next_id, 
                                           lambda n : (n['bigg_id'], n['name']), 
                                           x_off = j * block_width, y_off = i * block_height)

                        next_id = self.copy_map_reactions(to_draw.escher_map, 
                                                     master.escher_map, 
                                                     next_id, 
                                                     node_id_remap, compound_remap, 
                                                     lambda n : (n['bigg_id'], n['name']), 
                                                     x_off = j * block_width, y_off = i * block_height)
                    map_index += 1

        return master