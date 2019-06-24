class EscherMap:
    
    def __init__(self, escher_map):
        self.escher_map = escher_map
        self.escher_graph = escher_map[1]
        self.escher_data = escher_map[0]
        
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