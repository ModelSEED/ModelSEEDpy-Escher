import logging
logger = logging.getLogger(__name__)

def is_reversible(model_rxn):
    return model_rxn.lower_bound < 0 and model_rxn.upper_bound > 0
def reverse_stoichimetry(s):
    return dict(map(lambda x : (x[0], x[1] * -1), s.items()))

def validate_map_with_model(escher_model_map, model):
    model_spi_ids = set(map(lambda x : x.id, model.metabolites))
    model_rxn_ids = set(map(lambda x : x.id, model.reactions))
    for o in escher_model_map.metabolites:
        if o['bigg_id'] in model_spi_ids:
            model_spi = model.metabolites.get_by_id(o['bigg_id'])
            if not model_spi.name == o['name']:
                logger.debug("[species] (name) %s -> %s", o['name'], model_spi.name)
                o['name'] = model_spi.name
        else:
            logger.warning('[%s] not in map', o['bigg_id'])
    
    for map_rxn in escher_model_map.reactions:
        map_rxn_id = map_rxn['bigg_id']
        reversibility = map_rxn['reversibility']
        if map_rxn_id in model_rxn_ids:
            model_rxn = model.reactions.get_by_id(map_rxn_id)
            model_s = dict(map(lambda p : (p[0].id, p[1]), model_rxn.metabolites.items()))
            map_s = dict(map(lambda x : (x['bigg_id'], x['coefficient']), map_rxn['metabolites']))
            if model_s == map_s:
                if not reversibility == is_reversible(model_rxn):
                    map_rxn['reversibility'] = is_reversible(model_rxn)
                    logger.debug('%s incorrect reversiblity: %s, model: [%.2f, %.2f]', map_rxn_id, reversibility, model_rxn.lower_bound, model_rxn.upper_bound)
            elif reverse_stoichimetry(map_s) == model_s:
                logger.warning('%s contains LR/RL direction', map_rxn_id)
                for o in map_rxn['metabolites']:
                    o['coefficient'] = o['coefficient'] * - 1
                map_rxn['reversibility'] = is_reversible(model_rxn)
            else:
                matching_cpds = set(model_s) & set(map_s)
                print(map_rxn_id)
                print('map  :',map_s)
                print('model:', model_s)
        else:
            print(map_rxn_id)