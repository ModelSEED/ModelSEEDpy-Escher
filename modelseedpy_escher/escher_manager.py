from modelseedpy_escher.core.escher_map import EscherMap
import json

class EscherManager:
    
    def __init__(self, escher_instance):
        self.escher = escher_instance
        
    def list_datasets(self):
        res = set()
        for o in self.escher.list_cached_models():
            res.add(o['organism'])

        return res
    
    def list_models(self, dataset):
        res = set()
        for o in self.escher.list_cached_models():
            if o['organism'] == dataset:
                res.add(o['model_name'])

        return res

    def list_maps(self, dataset):
        res = set()
        for o in self.escher.list_cached_maps():
            if o['organism'] == dataset:
                res.add(o['map_name'])

        return res
    
    def get_map(self, dataset, model_id, map_id):
        escher_map = None
        with open("{}/maps/{}/{}.{}.json".format(self.escher.get_cache_dir(), dataset, model_id, map_id), 'r') as fh:
            escher_map = json.loads(fh.read())
        return EscherMap(escher_map)
    
    def save_map(self, dataset, model_id, map_id, escher_map):
        
        return False
    
    def save_model(self, dataset, model_id, escher_model):
        
        return False
