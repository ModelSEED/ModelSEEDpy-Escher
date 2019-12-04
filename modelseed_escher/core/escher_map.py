import json
import math
import copy
import escher

class EscherMap:
    
    def __init__(self, escher_map):
        self.escher_map = escher_map
        self.escher_graph = escher_map[1]
        self.escher_data = escher_map[0]
        
    def get_next_id(self):
        next_id = 0
        for node_id in self.escher_graph['nodes']:
            if int(node_id) >= next_id:
                next_id = int(node_id) + 1

        return next_id
        
    def swap_ids(self, cpd_remap, rxn_remap):
        for map_uid in self.escher_graph['nodes']:
            node = self.escher_graph['nodes'][map_uid]
            if node['node_type'] == 'metabolite' and node['bigg_id'] in cpd_remap:
                node['bigg_id'] = cpd_remap[node['bigg_id']]
        for map_uid in self.escher_graph['reactions']:
            map_reaction = self.escher_graph['reactions'][map_uid]
            if map_reaction['bigg_id'] in rxn_remap:
                map_reaction['bigg_id'] = rxn_remap[map_reaction['bigg_id']]
            for m in map_reaction['metabolites']:
                if m['bigg_id'] in cpd_remap:
                    m['bigg_id'] = cpd_remap[m['bigg_id']]
    
    def add_uid_to_reaction_metabolites(self):
        node_uid_map = {}
        for node_uid in self.escher_graph['nodes']:
            node = self.escher_graph['nodes'][node_uid]
            if node['node_type'] == 'metabolite':
                bigg_id = node['bigg_id']
                if not bigg_id in node_uid_map:
                    node_uid_map[bigg_id] = set()
                node_uid_map[bigg_id].add(node_uid)

        for node_uid in self.escher_graph['reactions']:
            rnode = self.escher_graph['reactions'][node_uid]
            met_id_to_uid = {}
            for s_uid in rnode['segments']:
                s = rnode['segments'][s_uid]
                from_node_id = rnode['segments'][s_uid]['from_node_id']
                to_node_id = rnode['segments'][s_uid]['to_node_id']
                from_node = self.escher_graph['nodes'][from_node_id]
                to_node = self.escher_graph['nodes'][to_node_id]
                if from_node['node_type'] == 'metabolite':
                    if not from_node['bigg_id'] in met_id_to_uid:
                        met_id_to_uid[from_node['bigg_id']] = from_node_id
                    else:
                        print('!!!', from_node['bigg_id'])
                if to_node['node_type'] == 'metabolite':
                    if not to_node['bigg_id'] in met_id_to_uid:
                        met_id_to_uid[to_node['bigg_id']] = to_node_id
                    else:
                        print('!!!', to_node['bigg_id'])

            for m in rnode['metabolites']:
                if m['bigg_id'] in met_id_to_uid:
                    m['node_uid'] = met_id_to_uid[m['bigg_id']]

            #print(met_id_to_uid)
            #print(rnode)
            #break

        return node_uid_map
    
    def delete_metabolites(self, cpd_ids):
        delete_uids = set()
        for map_uid in self.escher_graph['nodes']:
            node = self.escher_graph['nodes'][map_uid]
            if node['node_type'] == 'metabolite':
                node_id = node['bigg_id']

                if node_id in cpd_ids:
                    delete_uids.add(map_uid)

        for map_uid in delete_uids:
            del self.escher_graph['nodes'][map_uid]
            
    def delete_reactions(self, rxn_ids, remove_compounds = False):
        updated = {}
        delete_markers = set()
        delete_compounds = set()
        tagged_compounds = set()
        for map_uid in self.escher_graph['reactions']:
            rnode = self.escher_graph['reactions'][map_uid]
            if not rnode['bigg_id'] in rxn_ids:
                updated[map_uid] = rnode
            else:
                for s_uid in rnode['segments']:
                    s = rnode['segments'][s_uid]
                    from_node_id = s['from_node_id']
                    to_node_id = s['to_node_id']
                    n_type = self.escher_graph['nodes'][from_node_id]['node_type']
                    if n_type == 'metabolite':
                        delete_compounds.add(from_node_id)
                    else:
                        delete_markers.add(from_node_id)
                    #print(from_node_id, n_type)
                    n_type = self.escher_graph['nodes'][to_node_id]['node_type']
                    if n_type == 'metabolite':
                        delete_compounds.add(to_node_id)
                    else:
                        delete_markers.add(to_node_id)
                    #print(to_node_id, n_type)
                    #print(s_uid, s)
                #print(rnode)
        #delete also midmarkers
        for map_uid in delete_markers:
            del self.escher_graph['nodes'][map_uid]
        #delete oprhan compounds
        if remove_compounds:
            pass
        self.escher_graph['reactions'] = updated
    
    @property
    def nodes(self):
        return self.escher_graph['nodes']
    
    #@property
    #def reactions(self):
    #    return self.escher_graph['reactions']
    
    @property
    def reactions(self):
        reactions = list(
            self.escher_graph['reactions'].values()
        )

        return reactions
    
    @property
    def metabolites(self):
        metabolites = list(
            filter(
                lambda o : o['node_type'] == 'metabolite', 
                self.escher_graph['nodes'].values()
            )
        )

        return metabolites
    
    def generate_coords(self, func, n, x1, x2, x, y, orient, prod):
        a = x
        coords = []
        for pos_x in range(x1, x2, int((x2 - x1) / n)):
            print(pos_x + x, func(pos_x) + y)
            if prod:
                a -= 500
            coords.append([rotate(a, y, pos_x + x, func(pos_x) + y, orient)[0], rotate(a, y, pos_x + x, func(pos_x) + y, orient)[1]])
            if prod:
                a += 500
        return coords
    
    def get_nodes(self, stoich, rxn_numb):
        prim_node_ids = list()
        second_node_ids = list()
        node_ids2 = list()

        for k in stoich[str(rxn_numb)]:
            node_ids2.append(k)

        for k in range(0, len(node_ids2)):
            if stoich[str(rxn_numb)][node_ids2[k]] < 0:
                prim_node_ids.append(node_ids2.pop(k))
                break

        for k in range(0, len(node_ids2)):
            if stoich[str(rxn_numb)][node_ids2[k]] > 0:
                prim_node_ids.append(node_ids2.pop(k))
                break

        for k in node_ids2:
            second_node_ids.append(k)

        return prim_node_ids, second_node_ids
    
    def add_reaction(self, stoich, x_coord, y_coord, orient, reaction_length, prim_react, prim_prod, layout, rxn_numb, cc_coords, cc, prev_prod, prev_id):
        # base case
        if rxn_numb >= len(stoich):
            return
        # add midmarker node

        midmark = {'y': 0, 'x': 0, 'node_type': 'midmarker'}
        midmark_top = copy.deepcopy(midmark)
        midmark_bottom = copy.deepcopy(midmark)
        prim_react_id = 0
        prim_prod_id = 0 
        midmark_id = 0
        midmark_top_id = 0
        midmark_bottom_id = 0

        prim_nodes = self.get_nodes(stoich, rxn_numb)[0]
        second_nodes = self.get_nodes(stoich, rxn_numb)[1]
        print(prim_nodes)
        # find largest key in layout
        max = 0
        for k in layout['nodes'].keys():
            k = int(k)
            if k > max:
                max = k
        prim_react = {'node_is_primary': True,
                        'name': 'compound1',
                        'label_x': 0,
                        'node_type': 'metabolite',
                        'y': 0,
                        'x': 0,
                        'label_y': 0
                         }
        prim_prod = {'node_is_primary': True,
                        'name': 'compound1',
                        'label_x': 0,
                        'node_type': 'metabolite',
                        'y': 0,
                        'x': 0,
                        'label_y': 0
                         }  

        prim_react['bigg_id'] = prim_nodes[0]
        prim_prod['bigg_id'] = prim_nodes[1]

        # Check to see if either node is chained    
        if prim_react['bigg_id'] not in cc or (prev_prod == prim_react['bigg_id'] or prev_prod == None):
            # check to see if product is chained
            if prim_prod['bigg_id'] in cc and (prev_prod != None and prev_prod != prim_react['bigg_id']):
                orient += math.pi / 2
                orient = add_prim_nodes(prim_react, prim_prod, orient, cc_coords, reaction_length, midmark, layout)
                orient += math.pi
                layout['nodes'][str(max + 1)] = prim_react
                layout['nodes'][str(max + 2)] = midmark
                layout['nodes'][str(max + 3)] = midmark_top
                layout['nodes'][str(max + 4)] = midmark_bottom
                prim_prod_id = cc_coords[prim_prod['bigg_id']]['node_id']
                prim_react_id = str(max + 1)
                midmark_id = str(max + 2)
                midmark_top_id = str(max + 3)
                midmark_bottom_id = str(max + 4)
                max += 4
            # chained but linear    
            else:
                prim_react['x'] = x_coord
                prim_react['y'] = y_coord
                prim_react['label_x'] = x_coord - 30
                prim_react['label_y'] = y_coord - 30
                for n in layout['nodes']:
                    if layout['nodes'][n]['node_type'] == 'metabolite':
                        while (math.cos(orient) * reaction_length) + prim_react['x'] == layout['nodes'][n]['x'] and  - (math.sin(orient) * reaction_length) + prim_react['y'] == layout['nodes'][n]['y']:
                            orient += math.pi / 2
                prim_prod['x'] = (math.cos(orient) * reaction_length) + prim_react['x'] 
                prim_prod['y'] = - (math.sin(orient) * reaction_length) + prim_react['y'] 
                prim_prod['label_x'] = (math.cos(orient) * reaction_length) + prim_react['label_x']
                prim_prod['label_y'] = - (math.sin(orient) * reaction_length) + prim_react['label_y']
                midmark['x'] = (math.cos(orient) * (reaction_length / 2)) + prim_react['x'] 
                midmark['y'] = - (math.sin(orient) * (reaction_length / 2)) + prim_react['y']
                layout['nodes'][str(max + 1)] = midmark
                layout['nodes'][str(max + 2)] = prim_prod
                layout['nodes'][str(max + 3)] = midmark_top
                layout['nodes'][str(max + 4)] = midmark_bottom
                prim_prod_id = str(max + 2)
                prim_react_id = prev_id
                midmark_top_id = str(max + 3)
                midmark_bottom_id = str(max + 4)
                midmark_id = str(max + 1)
                max += 4

                if rxn_numb == 0:
                    max -= 4
                    layout['nodes'][str(max + 1)] = prim_react
                    layout['nodes'][str(max + 2)] = prim_prod
                    layout['nodes'][str(max + 3)] = midmark
                    layout['nodes'][str(max + 4)] = midmark_top
                    layout['nodes'][str(max + 5)] = midmark_bottom
                    prim_react_id = str(max + 1)
                    midmark_id = str(max + 3)
                    midmark_top_id = str(max + 4)
                    midmark_bottom_id = str(max + 5)
                    max += 5
        else:
            if prim_prod['bigg_id'] in cc_coords.keys():
                orient += math.pi / 2
                add_prim_nodes(prim_react, prim_prod, orient, cc_coords, reaction_length, midmark, layout)
                orient += math.pi
                layout['nodes'][str(max + 1)] = prim_react
                layout['nodes'][str(max + 2)] = midmark
                layout['nodes'][str(max + 3)] = midmark_top
                layout['nodes'][str(max + 4)] = midmark_bottom
                prim_prod_id = cc_coords[prim_prod['bigg_id']]['node_id']
                prim_react_id = str(max + 1)
                midmark_id = str(max + 2)
                midmark_top_id = str(max + 3)
                midmark_bottom_id = str(max + 4)
                max += 4

            else:
                orient += 3 * math.pi / 2
                add_prim_nodes(prim_prod, prim_react, orient, cc_coords, reaction_length, midmark, layout)
                layout['nodes'][str(max + 2)] = prim_prod
                layout['nodes'][str(max + 3)] = midmark
                layout['nodes'][str(max + 4)] = midmark_top
                layout['nodes'][str(max + 5)] = midmark_bottom
                prim_react_id = cc_coords[prim_react['bigg_id']]['node_id']
                prim_prod_id = str(max + 2)
                midmark_id = str(max + 3)
                midmark_top_id = str(max + 4)
                midmark_bottom_id = str(max + 5)
                max += 5



        # add coordinates of branched reactions
        if prim_prod['bigg_id'] in cc:
            cc_coords[prim_prod['bigg_id']]['x'] = prim_prod['x']
            cc_coords[prim_prod['bigg_id']]['y'] = prim_prod['y']
            cc_coords[prim_prod['bigg_id']]['node_id'] = prim_prod_id

        if prim_react['bigg_id'] in cc:
            cc_coords[prim_react['bigg_id']]['x'] = prim_react['x']
            cc_coords[prim_react['bigg_id']]['y'] = prim_react['y']
            cc_coords[prim_react['bigg_id']]['node_id'] = prim_react_id


        if rxn_numb == 0 and prim_react['bigg_id'] in cc:
            cc_coords[prim_react['bigg_id']]['x'] = prim_react['x']
            cc_coords[prim_react['bigg_id']]['y'] = prim_react['y']
            cc_coords[prim_react['bigg_id']]['node_id'] = prim_react_id

        # get largest key for reactions
        maxR = 0
        for k in layout['reactions'].keys():
            k = int(k)
            if k > maxR:
                maxR = k

        # get largest key for segments
        maxS = 0
        for k in layout['reactions'].keys():
            for j in layout['reactions'][k]['segments'].keys():
                if int(j) > maxS:
                    maxS = int(j)

        layout['reactions'][str(maxR + 1)] = {'name': str(rxn_numb),
              'bigg_id': str(rxn_numb),
              'segments': {},
              'genes': [],
              'reversibility': False,
              'metabolites': [],
              'label_x': midmark['x'] + 50 * math.sin(orient),
              'label_y': midmark['y'] + 50 * math.cos(orient),
              'gene_reaction_rule': ''}


        midmark_top['x'] = midmark['x'] - 25 * (math.cos(orient))
        midmark_bottom['x'] = midmark['x'] + 25 * (math.cos(orient))
        midmark_top['y'] = midmark['y'] + 25 * (math.sin(orient))
        midmark_bottom['y'] = midmark['y'] - 25 * (math.sin(orient))


        # alternate secondary node positions    
        ops = [operator.add, operator.sub]

        top = False
        p_to_s = 100
        mod = 0
        deg = 0
        # add secondary nodes to layout
        if len(second_nodes) > 0:
            maxS = add_secondary_nodes(second_nodes, orient, prim_react['x'], prim_react['y'], layout, stoich[str(rxn_numb)], midmark_top_id, midmark_bottom_id, midmark_id, midmark_top, midmark_bottom, midmark, maxR)
        maxS += 1
        #     for n in range(2, len(node_ids)):
    #         add_length = reaction_length
    #         if mod % 2 == 0:
    #             p_to_s += 60
    #             deg += math.pi / 30

    #         secondary = {'node_is_primary': False,
    #                       'name': 'compound1',
    #                       'label_x': 0,
    #                       'node_type': 'metabolite',
    #                       'y': 0,
    #                       'x': 0,
    #                       'bigg_id': node_ids[n],
    #                       'label_y': 0
    #                      }

    #         # secondary nodes are reactants
    #         react = False
    #         if int(stoich[str(rxn_numb)][str(node_ids[n])]) < 0:
    #             add_length = 0
    #             react = True

    #         # transformation matrix
    #         secondary['x'] = (ops[not(top)](abs(add_length - p_to_s * math.cos(math.pi / 6)) * (math.cos(orient)), (p_to_s * math.sin(math.pi / 6) * (math.sin(orient)))) + prim_react['x']) 
    #         secondary['y'] = (- ops[top](abs(add_length - p_to_s * math.cos(math.pi / 6)) * (math.sin(orient)), (p_to_s * math.sin(math.pi / 6) * (math.cos(orient)))) + prim_react['y']) 
    #         secondary['label_x'] = ops[top](secondary['x'] ,(30 * math.sin(orient)))
    #         secondary['label_y'] = ops[not(top)](secondary['y'] ,(30 * math.cos(orient)))
    #         layout['nodes'][str(max + 1)] = secondary
    #         add_segment(layout, midmark_top_id, str(max + 1), maxR, maxS, midmark, midmark_top, secondary, react, orient, n)
    #         layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': stoich[str(rxn_numb)][secondary['bigg_id']], 'bigg_id': secondary['bigg_id'] })
    #         maxS += 1
    #         max += 1
    #         top = not(top)
    #         mod += 1



        # add segments
        add_segment(layout, midmark_top_id, prim_react_id, maxR, maxS, midmark, midmark_top, prim_react, 1, orient, 0)
        layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': stoich[str(rxn_numb)][prim_react['bigg_id']], 'bigg_id': prim_react['bigg_id'] })
        add_segment(layout, prim_prod_id, midmark_bottom_id, maxR, maxS + 1, midmark, prim_prod, midmark_bottom, 0, orient, 0)
        layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': stoich[str(rxn_numb)][prim_prod['bigg_id']], 'bigg_id': prim_prod['bigg_id'] })
        add_segment(layout, midmark_id, midmark_top_id, maxR, maxS + 2, midmark, midmark, midmark_top, None, orient, 0)
        add_segment(layout, midmark_bottom_id, midmark_id, maxR, maxS + 3, midmark, midmark_bottom, midmark, None, orient, 0)


        # recursively place next reaction
        add_reaction(stoich, prim_prod['x'], prim_prod['y'], 3 * math.pi / 2, 500, 't1', 't2', layout, rxn_numb + 1, cc_coords, cc, prim_prod['bigg_id'], prim_prod_id )
    
    def rotate(self, x, y, x1, y1, orient):
        qx = x + math.cos(orient) * (x1 - x) - math.sin(orient) * (y1 - y)
        qy = y - (math.sin(orient) * (x1 - x) + math.cos(orient) * (y1 - y))
        return qx, qy
    
    def add_prim_nodes(self, prim_react, prim_prod, orient, cc_coords, reaction_length, midmark, layout):
        prim_prod['x'] = cc_coords[prim_prod['bigg_id']]['x']
        prim_prod['y'] = cc_coords[prim_prod['bigg_id']]['y']
        prim_prod['label_x'] = prim_prod['x'] - 30
        prim_prod['label_y'] = prim_prod['y'] - 30
        for n in layout['nodes']:
            if layout['nodes'][n]['node_type'] == 'metabolite':
                while (math.cos(orient) * reaction_length) + prim_prod['x'] == layout['nodes'][n]['x'] and  - (math.sin(orient) * reaction_length) + prim_prod['y'] == layout['nodes'][n]['y']:
                    orient += math.pi / 2
        prim_react['x'] = (math.cos(orient) * reaction_length) + prim_prod['x'] 
        prim_react['y'] = - (math.sin(orient) * reaction_length) + prim_prod['y'] 
        prim_react['label_x'] = (math.cos(orient) * reaction_length) + prim_prod['label_x']
        prim_react['label_y'] = - (math.sin(orient) * reaction_length) + prim_prod['label_y']
        midmark['x'] = prim_react['x'] - (math.cos(orient) * (reaction_length / 2))   
        midmark['y'] = prim_react['y'] + (math.sin(orient) * (reaction_length / 2)) 
        return orient
    
    def add_segment(self, layout, to_node_id, from_node_id, maxR, maxS, midmark, to_node, from_node, top, orient, n):
        dist1 = 50
        dist2 = 150
        if to_node['node_type'] == 'metabolite':
            if not(to_node['node_is_primary']):
                dist1 = 50 - (2 * n)
                dist2 = 150 - (2 * n)
        if from_node['node_type'] == 'metabolite':
            if not(from_node['node_is_primary']):
                dist1 = 50 - (2 * n)
                dist2 = 150 - (2 * n)
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)] = {}
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['to_node_id'] = to_node_id
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['from_node_id'] = from_node_id
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2'] = {}
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1'] = {}
        if to_node['node_type'] == 'midmarker' and from_node['node_type'] == 'midmarker':
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2'] = None
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1'] = None
        else:
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = midmark['x']
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = midmark['y']
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = midmark['y']
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = midmark['x']
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = ops[not(top)](midmark['y'], dist1 * math.sin(orient))
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = ops[not(top)](midmark['y'], dist2 * math.sin(orient))
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = ops[top](midmark['x'], dist1 * math.cos(orient))
            layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = ops[top](midmark['x'], dist2 * math.cos(orient))
    
    def add_secondary_nodes(self, second_nodes, orient, x_coord, y_coord, layout, stoich, midmark_top_id, midmark_bottom_id, midmark_id, midmark_top, midmark_bottom, midmark, maxR):
        react_coords = []
        react_nodes = []
        prod_coords = []
        for n in range(0, len(second_nodes) - 1):
            if stoich[second_nodes[n]] < 0:
                react_nodes.append(second_nodes.pop(n))

        prod_nodes = second_nodes 
        max = 0
        for n in layout['nodes']:
            if int(n) > max:
                max = int(n)
        maxS = 0
        for k in layout['reactions']:
            for j in layout['reactions'][k]['segments']:
                if int(j) > maxS:
                    maxS = int(j)

        if len(react_nodes) > 0: 
            react_coords = generate_coords(lambda x : pow(x, 2) / 30, len(react_nodes), 60, 160, x_coord, y_coord, orient, False)
        for c in range(0, len(react_coords)):
            secondary = {'node_is_primary': False,
                          'name': 'compound1',
                          'label_x': 0,
                          'node_type': 'metabolite',
                          'y': 0,
                          'x': 0,
                          'bigg_id': '',
                          'label_y': 0
                         }
            secondary['x'] = react_coords[c][0]
            secondary['label_x'] = react_coords[c][0]
            secondary['y'] = react_coords[c][1]
            secondary['label_y'] = react_coords[c][1] + 30
            secondary['bigg_id'] = react_nodes[c]
            layout['nodes'][str(max + 1)] = secondary
            add_segment(layout, str(max + 1), midmark_top_id, maxR, maxS, midmark, secondary, midmark_top, True, orient, c)
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': react_nodes[c]})
            max += 1
            maxS += 1

        if len(react_nodes) > 0:         
            prod_coords = generate_coords(lambda x : pow(x, 2) / 30, len(prod_nodes), -60, -160, x_coord + 500, y_coord, orient, True)
        for c in range(0, len(prod_coords) - 1):
            secondary = {'node_is_primary': False,
                          'name': 'compound1',
                          'label_x': 0,
                          'node_type': 'metabolite',
                          'y': 0,
                          'x': 0,
                          'bigg_id': '',
                          'label_y': 0
                         }
            secondary['x'] = prod_coords[c][0]
            secondary['label_x'] = prod_coords[c][0]
            secondary['y'] = prod_coords[c][1]
            secondary['label_y'] = prod_coords[c][1] + 30
            secondary['bigg_id'] = prod_nodes[c]
            layout['nodes'][str(max + 1)] = secondary
            add_segment(layout, str(max + 1), midmark_bottom_id, maxR, maxS, midmark, secondary, midmark_bottom, False, orient, c)
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': prod_nodes[c]})
            max += 1
            maxS += 1

        return maxS
    
    def add_curved_segment(self, layout, to_node_id, from_node_id, x, y):
        maxS = 0
        for k in layout['reactions']:
            for j in layout['reactions'][k]['segments']:
                if int(j) > maxS:
                    maxS = int(j)
        maxR = 0
        for k in layout['reactions'].keys():
            k = int(k)
            if k > maxR:
                maxR = k
        layout['reactions'][str(maxR + 1)] = {'name': 1,
                  'bigg_id': 1,
                  'segments': {},
                  'genes': [],
                  'reversibility': False,
                  'metabolites': [],
                  'label_x': 0,
                  'label_y': 0,
                  'gene_reaction_rule': ''}
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)] = {}
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['to_node_id'] = to_node_id
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['from_node_id'] = from_node_id
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2'] = {}
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1'] = {}
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = x - 30
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = y - 30
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = y - 30
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = x - 30

    def add_branches(self, layout, prim_react, prim_react_id, prim_nodes, numb, nodes_in_layout, reaction_length, max):
        maxS = 0
        for k in layout['reactions']:
            for j in layout['reactions'][k]['segments']:
                if int(j) > maxS:
                    maxS = int(j)
        maxR = 0
        for k in layout['reactions'].keys():
            k = int(k)
            if k > maxR:
                maxR = k
        x = 0
        midmark = {'y': prim_react['y'], 'x': prim_react['x'] + (reaction_length / 2), 'node_type': 'midmarker'}
        midmark_top = copy.deepcopy(midmark)
        midmark_bottom = copy.deepcopy(midmark)
        for n in prim_nodes.values():
            if n[0][0] == prim_react['bigg_id'] and n[1][0] not in nodes_in_layout:

                prim_prod = {'node_is_primary': True,
                        'name': 'compound1',
                        'label_x': prim_react['x'] + reaction_length,
                        'node_type': 'metabolite',
                        'y': prim_react['y'] - x * reaction_length,
                        'x': prim_react['x'] + reaction_length,
                        'label_y': prim_react['y'] - x * reaction_length,
                        'bigg_id': n[1][0]
                         }
                orient = math.atan((x * reaction_length) / (reaction_length))
                midmark = {'y': prim_prod['y'] + reaction_length / 2, 'x': (prim_prod['x'] - reaction_length / 2) - 100, 'node_type': 'midmarker'}
                midmark_top = {'y': prim_prod['y'] + reaction_length / 2 - (30 * math.sin(orient)), 'x': (prim_prod['x'] - reaction_length / 2 + (30 * math.cos(orient)) - 100) , 'node_type': 'midmarker'}
                midmark_bottom = {'y': prim_prod['y'] + reaction_length / 2 + (30 * math.sin(orient)), 'x': (prim_prod['x'] - reaction_length / 2 - (30 * math.cos(orient)) - 100), 'node_type': 'midmarker'}
                layout['nodes'][str(max + 1)] = prim_prod
                layout['nodes'][str(max + 2)] = midmark
                layout['nodes'][str(max + 3)] = midmark_top 
                layout['nodes'][str(max + 4)] = midmark_bottom
                add_curved_segment(layout, max + 1, max + 3, ((prim_prod['x'] - reaction_length / 2 + 10) + (prim_react['x'] + reaction_length)) / 2, ((prim_react['y'] - x * reaction_length) + (prim_react['y'] - x * reaction_length)) / 2 )
                add_curved_segment(layout, max + 4, prim_react_id, ((prim_prod['x'] - reaction_length / 2 - 10) + (prim_react['x'])) / 2, ((prim_prod['y'] + reaction_length / 2 + 10) + (prim_react['y'])) / 2 )
                add_segment(layout, max + 2, max + 4, maxR, maxS + 2, midmark, midmark, midmark_bottom, True, orient, 1)
                add_segment(layout, max + 3, max + 2, maxR, maxS + 3, midmark, midmark_top, midmark, True, orient, 1)
                maxS += 3
                max += 4
                x += 1
            if n[1][0] == prim_react['bigg_id'] and n[0][0] not in nodes_in_layout:
                prim_prod = {'node_is_primary': True,
                        'name': 'compound1',
                        'label_x': prim_react['x'] + reaction_length,
                        'node_type': 'metabolite',
                        'y': prim_react['y'] - x * reaction_length,
                        'x': prim_react['x'] + reaction_length,
                        'label_y': prim_react['y'] - x * reaction_length,
                        'bigg_id': n[0][0]
                         }
                orient = math.atan((x * reaction_length) / (reaction_length))
                midmark = {'y': prim_prod['y'] + reaction_length / 2, 'x': prim_prod['x'] - reaction_length / 2, 'node_type': 'midmarker'}
                midmark_top = {'y': prim_prod['y'] + reaction_length / 2 - (30 * math.sin(orient)), 'x': prim_prod['x'] - reaction_length / 2 + (30 * math.cos(orient)) , 'node_type': 'midmarker'}
                midmark_bottom = {'y': prim_prod['y'] + reaction_length / 2 + (30 * math.sin(orient)), 'x': prim_prod['x'] - reaction_length / 2 - (30 * math.cos(orient)), 'node_type': 'midmarker'}
                layout['nodes'][str(max + 1)] = prim_prod
                layout['nodes'][str(max + 2)] = midmark
                layout['nodes'][str(max + 3)] = midmark_top 
                layout['nodes'][str(max + 4)] = midmark_bottom
                add_curved_segment(layout, max + 1, max + 3, ((prim_prod['x'] - reaction_length / 2 + 10) + (prim_react['x'] + reaction_length)) / 2, ((prim_react['y'] - x * reaction_length) + (prim_react['y'] - x * reaction_length)) / 2 )
                add_curved_segment(layout, max + 4, prim_react_id, ((prim_prod['x'] - reaction_length / 2 - 10) + (prim_react['x'])) / 2, ((prim_prod['y'] + reaction_length / 2 + 10) + (prim_react['y'])) / 2 )
                add_segment(layout, max + 2, max + 4, maxR, maxS + 2, midmark, midmark, midmark_bottom, True, orient, 1)
                add_segment(layout, max + 3, max + 2, maxR, maxS + 3, midmark, midmark_top, midmark, True, orient, 1)
                maxS += 3
                max += 4
                x += 1
        return max
    
    def add_reactions(self, stoich, cc, x, y, layout):
        cc_coords = {}
        for c in cc:
            cc_coords[c] = {'x': 0, 'y': 0}
        self.add_reaction(stoich, x, y, 3 * math.pi / 2, 500, 't1', 't2', layout, 0, cc_coords, cc, None, None)
    
    def display_in_notebook(self):
        builder = escher.Builder(map_json=json.dumps(self.escher_map))
        return builder.display_in_notebook()
    
    def set_to_origin(self):
        offset_x = self.escher_graph['canvas']['x']
        offset_y = self.escher_graph['canvas']['y']
        self.escher_graph['canvas']['x'] = 0
        self.escher_graph['canvas']['y'] = 0

        for uid in self.nodes:
            n = self.nodes[uid]
            n['x'] -= offset_x
            n['y'] -= offset_y
            if 'label_x' in n:
                n['label_x'] -= offset_x
            if 'label_y' in n:
                n['label_y'] -= offset_y
        for uid in self.escher_graph['reactions']:
            rnode = self.escher_graph['reactions'][uid]
            rnode['label_x'] -= offset_x
            rnode['label_y'] -= offset_y
            for s_uid in rnode['segments']:
                segment = rnode['segments'][s_uid]
                if 'b1' in segment and segment['b1']:
                    segment['b1']['x'] -= offset_x
                    segment['b1']['y'] -= offset_y
                if 'b2' in segment and segment['b2']:
                    segment['b2']['x'] -= offset_x
                    segment['b2']['y'] -= offset_y
        for uid in self.escher_graph['text_labels']:
            tlabel = self.escher_graph['text_labels'][uid]
            tlabel['x'] -= offset_x
            tlabel['y'] -= offset_y
    
    def from_json(filename):
        escher_map = None
        with open(filename, 'r') as f:
            escher_map = EscherMap(json.loads(f.read()))
        return escher_map
      
    def clone(self):
        data = copy.deepcopy(self.escher_map)
        return EscherMap(data)