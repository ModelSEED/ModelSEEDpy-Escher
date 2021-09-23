import logging

logger = logging.getLogger(__name__)


def is_reversible(model_rxn):
    return model_rxn.lower_bound < 0 and model_rxn.upper_bound > 0


def reverse_stoichimetry(s):
    return dict(map(lambda x : (x[0], x[1] * -1), s.items()))


def validate_map_with_model(escher_model_map, model):
    model_spi_ids = set(map(lambda x: x.id, model.metabolites))
    model_rxn_ids = set(map(lambda x: x.id, model.reactions))
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
            logger.warning('[%s] not in map', map_rxn_id)


class RefitMap:

    def __init__(self, modelseed):
        self.ms = modelseed
        self.updated_node = set()

    def refit_map_metabolite(self, o, cmp_token):
        if o['bigg_id'] not in self.updated_node:
            o['bigg_id'] = "{}_{}".format(o['bigg_id'], cmp_token)
            o['name'] = "{} [{}]".format(o['name'], cmp_token)
            self.updated_node.add(o['bigg_id'])

    def refit_map_reaction(self, o, nodes, cmp_config):
        cmp_token = cmp_config['0']
        # print('refit_map_reaction', o['bigg_id'], cmp_config, cmp_token, len(cmp_token))
        if len(cmp_token) > 0:
            for s_uid in o['segments']:
                segment = o['segments'][s_uid]
                from_node = nodes[segment['from_node_id']]
                to_node = nodes[segment['to_node_id']]
                if to_node['node_type'] == 'metabolite':
                    self.refit_map_metabolite(to_node, cmp_token)
                if from_node['node_type'] == 'metabolite':
                    self.refit_map_metabolite(from_node, cmp_token)

            cmp_func = lambda x: {'bigg_id': x['bigg_id'] + '_' + cmp_token, 'coefficient': -1 * x['coefficient']}
            o['metabolites'] = list(map(cmp_func, o['metabolites']))

            o['bigg_id'] = "{}_{}".format(o['bigg_id'], cmp_token)
            o['name'] = "{} [{}]".format(o['name'], cmp_token)

    @staticmethod
    def reverse(map_reaction):
        rev_func = lambda x: {'bigg_id': x['bigg_id'], 'coefficient': -1 * x['coefficient']}
        metabolites = list(map(rev_func, map_reaction['metabolites']))
        map_reaction['metabolites'] = metabolites
        return map_reaction

    def refit(self, em, cmp_config):
        em = em.clone()
        self.updated_node = set()
        nodes = dict(map(lambda x: (x['uid'], x), em.nodes))
        seed_ids = {}
        for o in em.reactions:
            seed_id = o['bigg_id']
            rxn = self.ms.get_seed_reaction(seed_id)
            if rxn is not None:
                seed_ids[rxn.id] = rxn
                if rxn.data['direction'] == '=':
                    o['reversibility'] = True

                map_stoich = dict(map(lambda x: ((x['bigg_id'], '0'), x['coefficient']), o['metabolites']))
                cstoichiometry = rxn.cstoichiometry
                if rxn.data['direction'] == '<':
                    cstoichiometry = dict(map(lambda x: (x[0], -1 * x[1]), cstoichiometry.items()))
                if not map_stoich == cstoichiometry:
                    self.reverse(o)
                self.refit_map_reaction(o, nodes, cmp_config)
                if 'annotation' not in o:
                    o['annotation'] = {}
                if 'seed.reaction' not in o['annotation']:
                    o['annotation']['seed.reaction'] = []
                if type(o['annotation']['seed.reaction']) == str:
                    o['annotation']['seed.reaction'] = [o['annotation']['seed.reaction']]
                if rxn.id not in o['annotation']['seed.reaction']:
                    o['annotation']['seed.reaction'].append(rxn.id)
                o['annotation']['seed.compartment'] = ';'.join(map(lambda x: x[0] + ':' + x[1], cmp_config.items()))
        return em
