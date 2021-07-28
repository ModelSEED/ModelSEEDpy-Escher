def get_seed_ids(rxn):
    seed_ids = set()
    if '@' in rxn['id']:
        rxn_id, database = rxn['id'].split('@')
        seed_ids.add(rxn_id)
        if 'obsolete_seeds' in rxn:
            for rxn_id in rxn['obsolete_seeds']:
                rxn_id, database = rxn_id.split('@')
                seed_ids.add(rxn_id)
    else:
        if 'seed.reaction' in rxn['dblinks']:
            for rxn_id in rxn['dblinks']['seed.reaction']:
                seed_ids.add(rxn_id)
        if 'seed.obsolete' in rxn['dblinks']:
            for rxn_id in rxn['dblinks']['seed.obsolete']:
                seed_ids.add(rxn_id)
    if 'annotation' in rxn:
        if 'seed.reaction' in rxn['annotation']:
            seed_ids |= set(rxn['annotation']['seed.reaction'])
    return seed_ids

def get_cpd_seed_ids(model_metabolite):
    seed_ids = set()
    if '@' in model_metabolite['id']:
        seed_id, database = model_metabolite['id'].split('@')
        if database == 'SEED':
            seed_ids.add(seed_id)
            if 'obsolete_seeds' in model_metabolite:
                for seed_id in model_metabolite['obsolete_seeds']:
                    seed_id, database = seed_id.split('@')
                    seed_ids.add(seed_id)
    else:
        if 'cluster' in model_metabolite:
            for value in model_metabolite['cluster']:
                if '@' in value:
                    seed_id, database = value.split('@')
                    if database == 'SEED':
                        seed_ids.add(seed_id)
                else:
                    print('?', value)
                    
    if 'annotation' in model_metabolite:
        if 'seed.compound' in model_metabolite['annotation']:
            seed_ids |= set(model_metabolite['annotation']['seed.compound'])
    return seed_ids

class CobraJsonModel:
    
    def __init__(self, escher_model):
        self.escher_model = escher_model
        self.update_index()
    
    def update_index(self):
        self.id_to_index = {}
        for i in range(len(self.escher_model['metabolites'])):
            metabolite = self.escher_model['metabolites'][i]
            self.id_to_index[metabolite['id']] = i
            
    def get_metabolite_by_id(self, id):
        if id in self.id_to_index:
            return self.escher_model['metabolites'][self.id_to_index[id]]
        return None
    
    def asdasd(self):
        name_to_ids = {}
        #id_dups = {}
        for m_o in escher_model_original['metabolites']:
            name_to_ids[m_o['name']] = []
            #id_dups[m_o['id']] = []
        for m_o in escher_model_original['metabolites']:
            name_to_ids[m_o['name']].append(m_o['id'])
            #id_dups[m_o['id']].append(m_o)
        for k in id_dups:
            if len(id_dups[k]) > 1:
                1
                #print(k, id_dups[k])
                
    def merge_compounds(self, u_alias, name, ids, mapper):
        ids_mapped = set()

        for database in ids:
            if database in mapper:
                for cpd_id in ids[database]:
                    ids_mapped.add(cpd_id + '@' + mapper[database])
            else:
                print('ignore:', database)

        #print(ids_mapped)

        to_merge = []
        databases = set()
        for m in self.escher_model['metabolites']:
            database = 'universal'
            if '@' in m['id']:
                cpd_id, database = m['id'].split('@')
            if u_alias == m['id']:
                print('alias already exists:', m)
                return
            if m['id'] in ids_mapped:
                to_merge.append(m)

        if len(to_merge) == 0:
            print('nothing to merge')
            return

        merge_compound = {
            'charge': 0,
            'compartment': 'default',
            'formula': to_merge[0]['formula'],
            'id': u_alias,
            'name': name,
            'notes': {},
            'cluster' : []
        }

        for m in to_merge:
            merge_compound['cluster'].append(m['id'])

        merge_compound

        #add merge_compound
        #delete to_merge
        metabolites = []
        for m in self.escher_model['metabolites']:
            if not m['id'] in ids_mapped:
                metabolites.append(m)
        self.escher_model['metabolites'] = metabolites
        metabolites.append(merge_compound)

        for r in self.escher_model['reactions']:
            swap = ids_mapped & set(r['metabolites'].keys())
            if not 'replaced_cpds' in r:
                r['replaced_cpds'] = {}
            for cpd_id in swap:
                v = r['metabolites'][cpd_id]
                del r['metabolites'][cpd_id]
                r['metabolites'][u_alias] = v
                r['replaced_cpds'][u_alias] = cpd_id
    
    def set_annotation(self, id, bios_entry, bios_database):
        for m in self.escher_model['metabolites']:
            if id == m['id']:
                if not 'bios_references' in m:
                    m['bios_references'] = {}
                m['bios_references'][bios_database] = bios_entry
    
    def detect_models(self):
        model_ids = set()
        for m in self.escher_model['metabolites']:
            id = m['id']
            if '@' in id:
                model_id = id.split('@')[1]
                model_ids.add(model_id)
            else:
                model_ids.add('none')
        return model_ids
    
    #map_to target id
    #search ids to change
    def merge_model_nodes(self, map_to, search):
        #for m in escher_model['metabolites']:
        #    if m['id'] in search:
        #        #m['id'] = map_to
        #        print(m['id'], m['name'], '->', map_to)
        for r in self.escher_model['reactions']:
            replace = {}
            changed = False
            for id in r['metabolites']:
                if id in search and not id == map_to:
                    replace[map_to] = r['metabolites'][id]
                    changed = True
                else:
                    replace[id] = r['metabolites'][id]
            if changed:
                #print(r['id'], replace)
                r['metabolites'] = replace
                
    def merge_metabolites(self, keep, ids):
        if not keep in ids:
            print(keep, 'not in', ids)
            return None
        primary = self.get_metabolite_by_id(keep)
        if not 'bios_models' in primary:
            primary['bios_models'] = []
        primary['bios_models'].extend(ids)
        primary['bios_models'] = list(set(primary['bios_models']))
        to_delete = set()
        for id in ids:
            if not id == keep and id in self.id_to_index:
                to_delete.add(id)
            else:
                print('keep/skip:', id)

        self.delete_metabolites(to_delete)
        self.merge_model_nodes(keep, ids)
        return primary
    
    def delete_metabolites(self, ids):
        metabolites = []
        for m in self.escher_model['metabolites']:
            if not m['id'] in ids:
                metabolites.append(m)
        print('deleted', len(self.escher_model['metabolites']) - len(metabolites))
        self.escher_model['metabolites'] = metabolites
        self.update_index()
                
    def get_metabolite_groups(self):
        groups = {}
        for m in self.escher_model['metabolites']:
            bios_id = m['bios_id']
            m_id = m['id']
            if not m_id in groups:
                groups[m_id] = set()
            groups[m_id].add(bios_id)
        for i in groups:
            if len(groups[i]) > 1:
                #print(i, groups[i])
                model_count = {}
                for bios_id in groups[i]:
                    model_id = bios_id.split('@')[1]
                    if not model_id in model_count:
                        model_count[model_id] = 0
                    model_count[model_id] += 1
                    m = self.escher_model['metabolites'][bios_id_to_node_index[bios_id]]
                    bios_references = None
                    if 'bios_references' in m:
                        bios_references = m['bios_references']
                    if bios_references == None:
                        #print(i, bios_id_to_node_index[bios_id], bios_id, m['name'], bios_references)
                        1
                for model_id in model_count:
                    if model_count[model_id] > 1:
                        print(i, groups[i])
                        break
        return groups
    
    def initialize_bios_ids(self):
        for m in self.escher_model['metabolites']:
            m['bios_id'] = m['id']
        for r in self.escher_model['reactions']:
            r['bios_id'] = r['id']
            
    def map_escher_model_data(self, cpd_map_to_model, rxn_map_to_model):
        seed_match = set(rxn_map_to_model)
        rxn_remap = {}
        for model_reaction in self.escher_model['reactions']:
            seed_ids = get_seed_ids(model_reaction)
            match = seed_match & seed_ids
            if len(match) == 1:
                rxn_remap[model_reaction['id']] = rxn_map_to_model[match.pop()]
                #print(seed_ids, rxn_remap[model_reaction['id']])
            elif len(match) > 1:
                print('error', match, seed_ids)

        seed_match = set(cpd_map_to_model)
        cpd_remap = {}
        for model_metabolite in self.escher_model['metabolites']:
            seed_ids = get_cpd_seed_ids(model_metabolite)

            match = seed_match & seed_ids
            if len(match) == 1:
                cpd_remap[model_metabolite['id']] = cpd_map_to_model[match.pop()]
            elif len(match) > 1:
                cpd_remap[model_metabolite['id']] = cpd_map_to_model[sorted(match)[0]]
                print('error', match, seed_ids)
        return cpd_remap, rxn_remap