import logging
import json
import cobra
import escher
import cobrakbase
from modelseed_escher.core import CobraJsonModel
from modelseed_escher.core import EscherMap

logger = logging.getLogger(__name__)

def read_json(filename):
    data = None
    with open(filename, 'r') as f:
        data = json.loads(f.read())
    return data

def write_json(data, filename, pretty=False):
    with open(filename, 'w') as f:
        if pretty:
            f.write(json.dumps(data, indent=4, sort_keys=True))
        else:
            f.write(json.dumps(data))
            


class KBaseEscherViewer():
    
    def __init__(self, fbamodel, model_base, escher_model_file, fba = None, ESCHER_HOME = None, model_json = None):
        self.fbamodel = fbamodel
        self.fba = fba
        self.escher_base = model_base
        self.escher_model_file = escher_model_file
        if ESCHER_HOME == None:
            self.ESCHER_HOME = '/Users/fliu/workspace/jupyter/data/escher'
        else:
            self.ESCHER_HOME = ESCHER_HOME
        self.compartment_layer = {}
        self.map_list = list(filter(lambda d : d['organism'] == self.escher_base, escher.list_cached_maps()))
        
        self.maps = {}
        if not model_json == None:
            self.model_json = model_json
        else:
            cobra_model = cobrakbase.core.converters.KBaseFBAModelToCobraBuilder(self.fbamodel).build()
            self.model_json = cobra.io.to_json(cobra_model)
        
        self.flux_data = None
        if not self.fba == None:
            self.flux_data = self.load_flux(self.fba)
        
    def load_flux(self, fba):
        flux_data = {}
        for rxn_var in fba.data['FBAReactionVariables']:
            flux = rxn_var['value']
            rxn_id = rxn_var['modelreaction_ref'].split('/')[-1]
            flux_data[rxn_id] = flux
        for rxn_var in fba.data['FBABiomassVariables']:
            flux = rxn_var['value']
            rxn_id = rxn_var['biomass_ref'].split('/')[-1]
            flux_data[rxn_id + '_biomass'] = flux
        for rxn_var in fba.data['FBACompoundVariables']:
            flux = rxn_var['value']
            cpd_id = rxn_var['modelcompound_ref'].split('/')[-1]
            flux_data['EX_' + cpd_id] = -1 * flux
        return flux_data
        
    def load_model(self, compartment_match = ['c0']):
        escher_model_data = read_json(self.ESCHER_HOME + '/models/' + self.escher_base + '/' + self.escher_model_file)
        self.escher_model = CobraJsonModel(escher_model_data)
        
        for cmp_id in compartment_match:
            print(cmp_id)
            cpd_map_to_model, rxn_map_to_model = self.fetch_metabolites_and_reactions(self.fbamodel, cmp_id)
            cpd_remap, rxn_remap = self.escher_model.map_escher_model_data(cpd_map_to_model, rxn_map_to_model)

            self.compartment_layer[cmp_id] = {
                'cpd_map_to_model' : cpd_map_to_model,
                'rxn_map_to_model' : rxn_map_to_model,
                'cpd_remap' : cpd_remap,
                'rxn_remap' : rxn_remap,
            }
    
    @staticmethod
    def fetch_metabolites_and_reactions(fbamodel, compartment_match = 'c0'):
        rxn_map_to_model = {}

        for modelreaction in fbamodel.reactions:
            reaction_ref = modelreaction.data['reaction_ref'].split('/')[-1]
            seed_id = reaction_ref.split('_')[0]
            compartments = modelreaction.compartments
            logger.debug('%s[%s]', seed_id, compartments)
            if len(compartments) == 1 and compartments.pop() == compartment_match:
                rxn_map_to_model[seed_id] = modelreaction.id
            #print(modelreaction.id, seed_id, compartment)

        #compartment_match = 'c0'
        cpd_map_to_model = {}
        for metabolite in fbamodel.metabolites:
            seed_id = metabolite.data['compound_ref'].split('/')[-1]
            if not seed_id == 'cpd00000':
                compartment = metabolite.compartment
                logger.debug('%s[%s]', seed_id, compartment)
                if compartment == compartment_match:
                    if not seed_id in cpd_map_to_model:
                        cpd_map_to_model[seed_id] = metabolite.id
                    else:
                        logger.warning('double map %s -> %s, %s', seed_id, cpd_map_to_model[seed_id], metabolite.id)

        return cpd_map_to_model, rxn_map_to_model
        
    def build_escher_map(self, escher_map, model_json):
        map_json_str = json.dumps(escher_map.escher_map)
        builder = escher.Builder(map_json=map_json_str, model_json=model_json, reaction_data=self.flux_data)
        builder.set_highlight_missing(True)
        builder.set_enable_tooltips(True)
        builder.set_show_gene_reaction_rules(True)
        builder.set_and_method_in_gene_reaction_rule(True)
        #builder.set_hide_secondary_metabolites(True)
        return builder
    
    def generate_maps(self, layer):
        cpd_remap = self.compartment_layer[layer]['cpd_remap']
        rxn_remap = self.compartment_layer[layer]['rxn_remap']
        
        for o in self.map_list:
            #print(o)
            map_id = o['map_name']
            escher_map_data = read_json(self.ESCHER_HOME + '/maps/' + self.escher_base + '/' + map_id + '.json')
            escher_map = EscherMap(escher_map_data)
            rxn_in_map = set(map(lambda o : o['bigg_id'], escher_map.reactions))
            cpd_in_map = set(map(lambda o : o['bigg_id'], escher_map.metabolites))
            #print(map_id, len(cpd_in_map), len(rxn_in_map), len(rxn_in_map & set(rxn_remap)))
            escher_map.swap_ids(cpd_remap, rxn_remap)
            builder = self.build_escher_map(escher_map, self.model_json)
            self.maps[map_id] = [
                {
                    'metabolites' : len(cpd_in_map),
                    'reactions' : len(rxn_in_map),
                    'model_reactions' : len(rxn_in_map & set(rxn_remap)),
                    'model_genes' : '?',
                },
                builder
            ]
            #builder.save_html('../bin/escher/maps/' + map_id, overwrite=True)
            
    def save_maps(self, path):
        catalog = {}
        
        for map_id in self.maps:
            catalog[map_id] = {
                "src" : "maps/" + map_id + ".html"
            }
            catalog[map_id].update(self.maps[map_id][0])
            self.maps[map_id][1].save_html(path + '/' + map_id, overwrite=True)
        return catalog