import logging
logger = logging.getLogger(__name__)

def add_prefix(d, prefix_key, prefix_value):
    res = {}
    for k in d:
        k_mod = k if prefix_key == None else prefix_key + k
        res[k_mod] = set()
        for o in d[k]:
            res[k_mod].add(prefix_value + o)
    return res

def add_namespace(d, ns_key, ns_value):
    res = {}
    for k in d:
        k_ns = k if ns_key == None else k + '@' + ns_key
        res[k_ns] = set()
        for o in d[k]:
            res[k_ns].add(o + '@' + ns_value)
    return res

def get_mapping_from(modelseed_local, target_db):
    bigg_to_seed_rxn = {}
    bigg_to_seed = {}
    for seed_id in modelseed_local.compound_aliases:
        if target_db in modelseed_local.compound_aliases[seed_id]:
            for bigg_id in modelseed_local.compound_aliases[seed_id][target_db]:
                if not bigg_id in bigg_to_seed:
                    bigg_to_seed[bigg_id] = set()
                bigg_to_seed[bigg_id].add(seed_id)
    for seed_id in modelseed_local.reaction_aliases:
        rxn = modelseed_local.get_seed_reaction(seed_id)
        rxn_id = modelseed_local.get_non_obsolete(rxn)

        if target_db in modelseed_local.reaction_aliases[seed_id]:
            for bigg_id in modelseed_local.reaction_aliases[seed_id][target_db]:
                if not bigg_id in bigg_to_seed_rxn:
                    bigg_to_seed_rxn[bigg_id] = set()
                bigg_to_seed_rxn[bigg_id].add(rxn_id)
    return bigg_to_seed, bigg_to_seed_rxn

def get_mapping_to(modelseed_local, target_db):
    cpd_mapping = {}
    rxn_mapping = {}
    
    for seed_id in modelseed_local.compound_aliases:
        if target_db in modelseed_local.compound_aliases[seed_id]:
            for target_id in modelseed_local.compound_aliases[seed_id][target_db]:
                if not seed_id in cpd_mapping:
                    cpd_mapping[seed_id] = set()
                cpd_mapping[seed_id].add(target_id)
    for seed_id in modelseed_local.reaction_aliases:
        if target_db in modelseed_local.reaction_aliases[seed_id]:
            for target_id in modelseed_local.reaction_aliases[seed_id][target_db]:
                if not seed_id in rxn_mapping:
                    rxn_mapping[seed_id] = set()
                rxn_mapping[seed_id].add(target_id)
    return cpd_mapping, rxn_mapping


    
def remap_map_compounds(em, bigg_to_seed, id_function = lambda x : (x, 0)):
    map_compound_remap = {}
    unmaped = set()
    node_uid_cmp = {}
    for map_uid in em.escher_graph['nodes']:
        node = em.escher_graph['nodes'][map_uid]
        if node['node_type'] == 'metabolite':
            node_id = node['bigg_id']
            bigg_id, cmp = id_function(node_id)
            #bigg_id = node_id[:-2]
            #cmp = node_id[-1:]
            
            #print(node_id, bigg_id, cmp)
            if bigg_id in bigg_to_seed:
                node_uid_cmp[map_uid] = cmp
                map_compound_remap[node_id] = set(bigg_to_seed[bigg_id])
                #print(map_uid, bigg_id, cmp, bigg_to_seed[bigg_id])
            else:
                unmaped.add(node_id)
    return map_compound_remap, unmaped, node_uid_cmp

ALLOW_MISSING = set()

def remap_map_reactions(em, bigg_to_seed_rxn, map_compound_remap, get_rxn, to_match_func = lambda x : x):
    unmaped_rxn = set()
    map_reaction_remap = {}
    for map_uid in em.escher_graph['reactions']:
        rnode = em.escher_graph['reactions'][map_uid]
        node_id = rnode['bigg_id']
        if node_id in bigg_to_seed_rxn:
            for db_id in bigg_to_seed_rxn[node_id]:
                map_stoichiometry = get_stoichiometry(rnode)
                #print(db_id)
                #NEED KEGG/METACYC/BIGG provenance
                rxn_cstoich = get_rxn(db_id)
                
                #print(rxn_cstoich)
                if rxn_cstoich != None:
                    to_match = set(map(lambda x : to_match_func(x[0]), rxn_cstoich.keys()))
                    
                    #print(to_match, map_stoichiometry)
                    mapping, missing, to_match = match2(to_match, map_stoichiometry, map_compound_remap)
                    missing -= ALLOW_MISSING
                    #print(missing)
                    if len(missing) == 0:
                        #print(map_uid, node_id, db_id, map_stoichiometry, rxn.cstoichiometry)
                        #print(mapping, missing, to_match)
                        if not node_id in map_reaction_remap:
                            map_reaction_remap[node_id] = set()
                        map_reaction_remap[node_id].add(db_id)
                        for s_id in rnode['segments']:
                            s = rnode['segments'][s_id]
                            #print(s_id, s)
                            break
                    else:
                        print(node_id, db_id, missing)
                        unmaped_rxn.add(node_id)
        else:
            print(node_id, "?")
            unmaped_rxn.add(node_id)
    return map_reaction_remap, unmaped_rxn

def match(rxn, map_stoichiometry, map_compound_remap):
    to_match = set(map(lambda x : x[0], rxn.cstoichiometry.keys()))
    missing = set()
    mapping = {}
    for cpd_id in map_stoichiometry:
        if cpd_id in map_compound_remap:
            for other_id in map_compound_remap[cpd_id]:
                if other_id in to_match:
                    mapping[cpd_id] = other_id
    missing =  to_match - set(mapping.values())
    return mapping, missing, to_match

def match2(to_match, map_stoichiometry, map_compound_remap):
    #to_match = set(map(lambda x : x[0], rxn.cstoichiometry.keys()))
    missing = set()
    mapping = {}
    for cpd_id in map_stoichiometry:
        if cpd_id in map_compound_remap:
            for other_id in map_compound_remap[cpd_id]:
                if other_id in to_match:
                    mapping[cpd_id] = other_id
    missing =  to_match - set(mapping.values())
    return mapping, missing, to_match

def get_stoichiometry(rnode):
    stoichiometry = {}
    for m in rnode['metabolites']:
        stoichiometry[m['bigg_id']] = m['coefficient']
        
    return stoichiometry

def add_compartment(em, node_uid_cmp):
    for node_uid in node_uid_cmp:
        if node_uid in em.escher_graph['nodes']:
            node = em.escher_graph['nodes'][node_uid]
            node['compartment'] = node_uid_cmp[node_uid]
            

def get_cstoich_list(em, rnode):
    cstoich_list = []
    
    for s_uid in rnode['segments']:
        s = rnode['segments'][s_uid]
        from_node_id = rnode['segments'][s_uid]['from_node_id']
        to_node_id = rnode['segments'][s_uid]['to_node_id']
        from_node = em.escher_graph['nodes'][from_node_id]
        to_node = em.escher_graph['nodes'][to_node_id]
        if from_node['node_type'] == 'metabolite':
            cstoich_list.append((from_node['bigg_id'], from_node['compartment']))
        if to_node['node_type'] == 'metabolite':
            cstoich_list.append((to_node['bigg_id'], to_node['compartment']))
            
    return cstoich_list