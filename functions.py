import math 
import copy
def add_secondary_nodes(second_react, second_prod, orient, x_coord, y_coord, layout, midmark_top_id, midmark_bottom_id, midmark_id, midmark_top, midmark_bottom, midmark, maxR, prim, reaction_length, rxn_numb):
    react_coords = []
    react_nodes = []
    prod_coords = []
    prod_nodes = []
    for n in second_react:
        react_nodes.append(n[0])
    for n in second_prod:
        prod_nodes.append(n[0])
            
    max = 0
    for n in layout['nodes']:
        if int(n) > max:
            max = int(n)    
    maxS = 0
    for k in layout['reactions']:
        for j in layout['reactions'][k]['segments']:
            if int(j) > maxS:
                maxS = int(j)
    
    if prim:
        temp = []
        temp = copy.deepcopy(react_nodes)
        react_nodes = copy.deepcopy(prod_nodes)
        prod_nodes = copy.deepcopy(temp)

        
        
    if len(react_nodes) > 0: 
        react_coords = generate_coords(lambda x : - pow(x, 2) / 50, len(react_nodes), 80, 160, x_coord, y_coord, 1)
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
        secondary['label_x'] = react_coords[c][0] - 90
        secondary['y'] = react_coords[c][1]
        secondary['label_y'] = react_coords[c][1] - 20       
        secondary['bigg_id'] = react_nodes[c]  
        layout['nodes'][str(max + 1)] = secondary
        add_segment(layout, str(max + 1), midmark_top_id, maxR, maxS, midmark, secondary, midmark_top, True, orient, c)
        if prim:
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': 1, 'bigg_id': react_nodes[c]})
        else:
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': react_nodes[c]})

        max += 1
        maxS += 1
      
    if len(prod_nodes) > 0:         
        prod_coords = generate_coords(lambda x : - pow(x, 2) / 50, len(prod_nodes), -60, -160, x_coord + reaction_length, y_coord, -1)
    for c in range(0, len(prod_coords)):
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
        secondary['label_x'] = prod_coords[c][0] - 25
        secondary['y'] = prod_coords[c][1]
        secondary['label_y'] = prod_coords[c][1] - 20
        secondary['bigg_id'] = prod_nodes[c]
        layout['nodes'][str(max + 1)] = secondary
        add_segment(layout, str(max + 1), midmark_bottom_id, maxR, maxS, midmark, secondary, midmark_bottom, False, orient, c)
        if prim:
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': prod_nodes[c]})
        else:
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': 1, 'bigg_id': prod_nodes[c]})
        max += 1
        maxS += 1
    
    return maxS, max
    

def add_curved_segment(layout, to_node_id, from_node_id, to_node, from_node, orient, maxR, x):
    maxS = 0
    for k in layout['reactions']:
        for j in layout['reactions'][k]['segments']:
            if int(j) > maxS:
                maxS = int(j)

    if x == 0:
        x = 1
    layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)] = {}
    layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['to_node_id'] = to_node_id
    layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['from_node_id'] = from_node_id
    layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2'] = {}
    layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1'] = {}
    if to_node['node_type'] == 'midmarker':
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = to_node['x'] - 100 * x
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = to_node['y']
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = from_node['y'] 
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = from_node['x'] + 100
    else:
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = to_node['x'] - 100 
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = to_node['y'] 
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = from_node['y']
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = from_node['x'] + 100
        
def get_max(layout):
    max = 0
    for n in layout['nodes']:
        if int(n) > max:
            max = int(n)
    return max
def generate_coords(func, n, x1, x2, x, y, cap):
    coords = []
    for pos_x in range(x1, x2 - (cap * (n)), int((x2 - x1) / n)):     
        coords.append([pos_x + x, func(pos_x) + y])
    return coords

def add_segment(layout, to_node_id, from_node_id, maxR, maxS, midmark, to_node, from_node, react, orient, n):
    import operator
    import math
    ops = [operator.add, operator.sub]
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
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = ops[not(react)](midmark['y'], dist1 * math.sin(orient))
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = ops[not(react)](midmark['y'], dist2 * math.sin(orient))
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = ops[react](midmark['x'], dist1 * math.cos(orient))
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = ops[react](midmark['x'], dist2 * math.cos(orient))
        
def add_segment_prim(layout, to_node_id, from_node_id, maxR, maxS, midmark, to_node, from_node, top, orient, n):
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
    if to_node['node_type'] == 'midmarker':
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = to_node['x'] - 200
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = to_node['y']
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = from_node['y'] - 100 * math.sin(orient) 
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = from_node['x'] + 100 * math.cos(orient) 
    else:
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['x'] = to_node['x'] - 100 
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b2']['y'] = to_node['y'] 
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['y'] = from_node['y']
        layout['reactions'][str(maxR + 1)]['segments'][str(maxS + 1)]['b1']['x'] = from_node['x'] + 100
        
    
def add_reaction(layout, prim_react, prim_react_id, prim_nodes, x, y, reaction_length, orient, rxn_numb, cc, nodes_in_layout, second_react_nodes, second_prod_nodes, prim_id):
    if len(prim_nodes) == 0:
        return
    max = 0
    for k in layout['nodes'].keys():
        k = int(k)
        if k > max:
            max = k
    if cc[prim_react['bigg_id']] > 2:
        max = add_branches(layout, prim_react, prim_react_id, prim_nodes, 4, nodes_in_layout, reaction_length, max, second_react_nodes, second_prod_nodes, rxn_numb, cc)
    else:
        prim_prod = {'node_is_primary': True,
                        'name': 'compound1',
                        'label_x': 0,
                        'node_type': 'metabolite',
                        'y': 0,
                        'x': 0,
                        'label_y': 0
                         }
        found = False
        coefficient = 0
        for key,reaction in prim_nodes.items():
            if reaction[0][0] == prim_react['bigg_id']:
                prim_prod['bigg_id'] = reaction[1][0]
                second_react = second_react_nodes[key]
                second_prod = second_prod_nodes[key]
                del prim_nodes[key]
                found = True
                coefficient = reaction[1][1]
                break
        
        if not(found):
            for key,reaction in prim_nodes.items():
                if reaction[1][0] == prim_react['bigg_id']:
                    prim_prod['bigg_id'] = reaction[0][0]
                    second_react = second_react_nodes[key]
                    second_prod = second_prod_nodes[key]
                    del prim_nodes[key]
                    found = True
                    coefficient = reaction[0][1]
                    break
                
        
        if not(found):
            return
        x = 0
        ort = 0
        for m in layout['nodes'].values():
                while m['x'] == prim_react['x'] + reaction_length and m['y'] == prim_react['y'] - reaction_length * x:
                    x += 1
                    ort = math.pi / 2
        prim_prod['x'] = prim_react['x'] + reaction_length
        prim_prod['y'] = prim_react['y'] - reaction_length * x
        prim_prod['label_x'] = prim_prod['x']
        prim_prod['label_y'] = prim_prod['y'] + 35

        
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
        midmark = {'y': prim_prod['y'], 'x': prim_prod['x'] - (reaction_length / 2), 'node_type': 'midmarker'}
        midmark_top = copy.deepcopy(midmark)
        midmark_bottom = copy.deepcopy(midmark)   
        midmark_top['x'] = midmark['x'] - 25 * (math.cos(orient))
        midmark_bottom['x'] = midmark['x'] + 25 * (math.cos(orient))
        midmark_top['y'] = midmark['y'] + 25 * (math.sin(orient))
        midmark_bottom['y'] = midmark['y'] - 25 * (math.sin(orient))
        if rxn_numb == 0:
            layout['nodes'][str(max + 1)] = prim_react
            prim_id = max + 1
        layout['nodes'][str(max + 2)] = prim_prod
        layout['nodes'][str(max + 3)] = midmark
        layout['nodes'][str(max + 4)] = midmark_top
        layout['nodes'][str(max + 5)] = midmark_bottom

        layout['reactions'][str(maxR + 1)] = {'name': str(rxn_numb),
              'bigg_id': str(rxn_numb),
              'segments': {},
              'genes': [],
              'reversibility': False,
              'metabolites': [],
              'label_x': midmark['x'] + 50 * math.sin(orient),
              'label_y': midmark['y'] + 50 * math.cos(orient),
              'gene_reaction_rule': ''}

        midmark_id = max + 3
        midmark_top_id = max + 4
        midmark_bottom_id = max + 5
        midmark_top['x'] = midmark['x'] - 25 * (math.cos(orient))
        midmark_bottom['x'] = midmark['x'] + 25 * (math.cos(orient))
        midmark_top['y'] = midmark['y'] + 25 * (math.sin(orient))
        midmark_bottom['y'] = midmark['y'] - 25 * (math.sin(orient))
        add_segment_prim(layout, max + 4, prim_id, maxR, maxS + 1, midmark, midmark_top, prim_react, 1, ort, 0)
        add_segment(layout, max + 3, max + 4, maxR, maxS + 2, midmark, midmark, midmark_top, 1, orient, 0)
        add_segment(layout, max + 5, max + 3, maxR, maxS + 3, midmark, midmark_bottom, midmark, 1, orient, 0)
        add_segment_prim(layout, max + 2, max + 5, maxR, maxS + 4, midmark, prim_prod, midmark_bottom, 1, orient, 0)
        prim = True
        if coefficient > 0:
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': 1, 'bigg_id': prim_prod['bigg_id']}) 
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': prim_react['bigg_id']})         
            prim = False        
        else:
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': prim_prod['bigg_id']}) 
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': 1, 'bigg_id': prim_react['bigg_id']}) 
       
        nodes_in_layout.append(prim_react)
        nodes_in_layout.append(prim_prod)
        add_secondary_nodes(second_react, second_prod, orient, prim_react['x'], prim_react['y'], layout, midmark_top_id, midmark_bottom_id, midmark_id, midmark_top, midmark_bottom, midmark, maxR, prim, reaction_length, rxn_numb)
        add_reaction(layout, prim_prod, max + 2, prim_nodes, x + reaction_length, y, reaction_length, orient, rxn_numb + 1, cc, nodes_in_layout, second_react_nodes, second_prod_nodes, max + 2)
def add_branches(layout, prim_react, prim_react_id, prim_nodes, numb, nodes_in_layout, reaction_length, max, second_react_nodes, second_prod_nodes, rxn_numb, cc):
    maxS = 0
    for k in layout['reactions']:
        for j in layout['reactions'][k]['segments']:
            if int(j) > maxS:
                maxS = int(j)
   
    x = 0
    for key, n in prim_nodes.items():
        if n[0][0] == prim_react['bigg_id'] and n[1][0] not in nodes_in_layout:
            max = get_max(layout)
            for m in layout['nodes'].values():
                while m['x'] == prim_react['x'] + reaction_length and m['y'] == prim_react['y'] - reaction_length * x:
                    x += 1                  
            prim_prod = {'node_is_primary': True,
                    'name': 'compound1',
                    'label_x': prim_react['x'] + reaction_length,
                    'node_type': 'metabolite',
                    'y': prim_react['y'] - reaction_length * x,
                    'x': prim_react['x'] + reaction_length,
                    'label_y': prim_react['y'] - reaction_length * x + 35,
                    'bigg_id': n[1][0]
                     }
            maxR = 0
            for k in layout['reactions'].keys():
                k = int(k)
                if k > maxR:
                    maxR = k
            orient = 0
            midmark = {'y': prim_react['y'] - x * reaction_length, 'x': prim_react['x'] + reaction_length / 2, 'node_type': 'midmarker'}
            midmark_top = {'y': prim_react['y'] - x * reaction_length, 'x': prim_react['x'] + reaction_length / 2 + 50, 'node_type': 'midmarker'}
            midmark_bottom = {'y': prim_react['y'] - x * reaction_length, 'x': prim_react['x'] + reaction_length / 2 - 50, 'node_type': 'midmarker'}
            layout['reactions'][str(maxR + 1)] = {'name': rxn_numb + 1,
              'bigg_id': rxn_numb,
              'segments': {},
              'genes': [],
              'reversibility': False,
              'metabolites': [],
              'label_x': midmark['x'] ,
              'label_y': midmark['y'] - 30,
              'gene_reaction_rule': ''}
            layout['nodes'][str(max + 1)] = prim_prod
            layout['nodes'][str(max + 2)] = midmark
            layout['nodes'][str(max + 3)] = midmark_top 
            layout['nodes'][str(max + 4)] = midmark_bottom
            prim_prod_id = max + 1
            midmark_id = max + 2
            midmark_top_id = max + 3
            midmark_bottom_id = max + 4
            second_react = second_react_nodes[key]
            second_prod = second_prod_nodes[key]
            add_curved_segment(layout, prim_prod_id, midmark_top_id, prim_prod, midmark_top, orient, maxR, x)
            add_curved_segment(layout, midmark_bottom_id, prim_react_id, midmark_bottom, prim_react, orient, maxR, x)
            maxes = add_secondary_nodes(second_react, second_prod, 0, midmark['x'] - reaction_length / 2, midmark['y'], layout, max + 3, max + 4, max + 2, midmark_top, midmark_bottom, midmark, maxR, False, reaction_length, rxn_numb)
            maxS = maxes[0]
            max = maxes[1]
            add_segment(layout, midmark_id, midmark_bottom_id, maxR, maxS + 1, midmark, midmark, midmark_bottom, True, orient, 1)
            add_segment(layout, midmark_top_id, midmark_id, maxR, maxS + 2, midmark, midmark_top, midmark, True, orient, 1)
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': 1, 'bigg_id': prim_prod['bigg_id']})                        
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': prim_react['bigg_id']}) 
            maxS += 2
            max += 4
            x += 1
            prim_nodes2 = copy.deepcopy(prim_nodes)
            del prim_nodes2[key]
            nodes_in_layout.append(prim_prod['bigg_id'])            
            add_reaction(layout, prim_prod, prim_prod_id, prim_nodes2, prim_prod['x'] + 500, prim_prod['y'], reaction_length, 0, rxn_numb + 1, cc, nodes_in_layout, second_react_nodes, second_prod_nodes, prim_prod_id)
            rxn_numb += 1
        if n[1][0] == prim_react['bigg_id'] and n[0][0] not in nodes_in_layout:
            max = get_max(layout)
            for m in layout['nodes'].values():
                while m['x'] == prim_react['x'] + reaction_length and m['y'] == prim_react['y'] - reaction_length * x:
                    x += 1 
            prim_prod = {'node_is_primary': True,
                    'name': 'compound1',
                    'label_x': prim_react['x'] + reaction_length,
                    'node_type': 'metabolite',
                    'y': prim_react['y'] - reaction_length * x,
                    'x': prim_react['x'] + reaction_length,
                    'label_y': prim_react['y'] - reaction_length * x + 35,
                    'bigg_id': n[0][0]
                     }
            maxR = 0
            for k in layout['reactions'].keys():
                k = int(k)
                if k > maxR:
                    maxR = k
            orient = 0
            midmark = {'y': prim_react['y'] - x * reaction_length, 'x': prim_react['x'] + reaction_length / 2, 'node_type': 'midmarker'}
            midmark_top = {'y': prim_react['y'] - x * reaction_length, 'x': prim_react['x'] + reaction_length / 2 + 50, 'node_type': 'midmarker'}
            midmark_bottom = {'y': prim_react['y'] - x * reaction_length, 'x': prim_react['x'] + reaction_length / 2 - 50, 'node_type': 'midmarker'}                      
            layout['reactions'][str(maxR + 1)] = {'name': rxn_numb + 1,
              'bigg_id': rxn_numb,
              'segments': {},
              'genes': [],
              'reversibility': False,
              'metabolites': [],
              'label_x': midmark['x'],
              'label_y': midmark['y'] - 30,
              'gene_reaction_rule': ''}
            layout['nodes'][str(max + 1)] = prim_prod
            layout['nodes'][str(max + 2)] = midmark
            layout['nodes'][str(max + 3)] = midmark_top 
            layout['nodes'][str(max + 4)] = midmark_bottom
            prim_prod_id = max + 1
            midmark_id = max + 2
            midmark_top_id = max + 3
            midmark_bottom_id = max + 4
            second_react = second_react_nodes[key]
            second_prod = second_prod_nodes[key]
            add_curved_segment(layout, prim_prod_id, midmark_top_id, prim_prod, midmark_top, orient, maxR, x)
            add_curved_segment(layout, midmark_bottom_id, prim_react_id, midmark_bottom, prim_react, orient, maxR, x)
            maxes = add_secondary_nodes(second_react, second_prod, 0, midmark['x'] - reaction_length / 2, midmark['y'], midmark['x'], midmark['y'], layout, stoich, max + 3, max + 4, max + 2, midmark_top, midmark_bottom, midmark, maxR, True, reaction_length, rxn_numb)
            maxS = maxes[0]
            max = maxes[1]
            add_segment(layout, midmark_id, midmark_bottom_id, maxR, maxS + 1, midmark, midmark, midmark_bottom, True, orient, 1)
            add_segment(layout, midmark_top_id, midmark_id, maxR, maxS + 2, midmark, midmark_top, midmark, True, orient, 1)
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': -1, 'bigg_id': prim_prod['bigg_id']})             
            layout['reactions'][str(maxR + 1)]['metabolites'].append({'coefficient': 1, 'bigg_id': prim_react['bigg_id']})             
            maxS += 2
            max += 4
            x += 1
            prim_nodes2 = copy.deepcopy(prim_nodes)
            del prim_nodes2[key]
            nodes_in_layout.append(prim_prod['bigg_id'])
            add_reaction(layout, prim_prod, prim_prod_id, prim_nodes2, prim_prod['x'] + 500, prim_prod['y'], reaction_length, 0, rxn_numb + 1, cc, nodes_in_layout, second_react_nodes, second_prod_nodes, prim_prod_id)
            rxn_numb += 1
    return max

def add_reactions(layout, stoich, x, y, reaction_length, orient, bigg_id):
    prim_react = {'node_is_primary': True,
                    'name': 'compound1',
                    'label_x': x,
                    'node_type': 'metabolite',
                    'y': y,
                    'x': x,
                    'label_y': y + 35,
                    'bigg_id': bigg_id
                     }
    prim_nodes = {}
    second_react_nodes = {}
    second_prod_nodes = {}
    for r in stoich:
        prim_nodes[r] = []
        second_react_nodes[r] = []
        second_prod_nodes[r] = []
        numb = 0
        for n in stoich[r]:
            if stoich[r][n] < 0 and numb == 0:
                dict = []
                dict.append(n)
                dict.append(stoich[r][n])
                prim_nodes[r].append(dict)
                numb +=1
            elif stoich[r][n] < 0:
                dict1 = []
                dict1.append(n)
                dict1.append(stoich[r][n])
                second_react_nodes[r].append(dict1)
        numb = 0
        for n in stoich[r]:
            if stoich[r][n] > 0 and numb == 0:
                dict2 = []
                dict2.append(n)
                dict2.append(stoich[r][n])
                prim_nodes[r].append(dict2)
                numb +=1
            elif stoich[r][n] > 0:
                dict3 = []
                dict3.append(n)
                dict3.append(stoich[r][n])
                second_prod_nodes[r].append(dict3)
    #count number of times compound appears in dict                
    nodes = []
    for r in prim_nodes.values():
        for n in r:
            if n[0] not in nodes:
                nodes.append(n[0])

    cc = {}
    for n in nodes:
        cc[n] = 0
    for n in nodes:
        for r in prim_nodes.values():
            for node in r[0]:
                if n == node:
                    cc[n] += 1
            for node in r[1]:
                if n == node:
                    cc[n] += 1


    nodes_in_layout = []
    add_reaction(layout, prim_react, 1, prim_nodes, x, y, reaction_length, orient, 0, cc, nodes_in_layout, second_react_nodes, second_prod_nodes, 0)
    
