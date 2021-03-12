import copy
import logging

logger = logging.getLogger(__name__)


def get_cmp_token(cmps):
    if 'k' in cmps:
        return 'k'
    if 'km' in cmps:
        return 'km'
    if len(cmps) == 1:
        return list(cmps)[0]
    if len(cmps) == 2:
        if 'b' in cmps and 'e' in cmps:
            return 'b'
        if 'e' in cmps and 'c' in cmps:
            return 'c'
        if 'e' in cmps and 'p' in cmps:
            return 'e'
        if 'c' in cmps:
            return list(filter(lambda x: not x == 'c', cmps))[0]

    return None


def modelseed_to_cobra_reaction(rxn, compartment=None, rxn_cmp_suffix='c',
                                lb=None, ub=None, rev='direction',
                                cs_override=None):
    if compartment is None:
        compartment = {'0': ''}
    crxn = {
        'id': 'RXN',
        'name': 'Reaction',
        'metabolites': {},
        'lower_bound': -1000.0,
        'upper_bound': 1000.0,
        'gene_reaction_rule': '',
        'annotation': {}
    }

    metabolites = {}
    cpd_cmp = {}
    # cs = get_cstoichiometry_plant(rxn.data)
    cs = rxn.cstoichiometry
    if not cs_override == None:
        logger.debug("cs_override %s -> %s", cs, cs_override)
        cs = cs_override

    logger.debug("%s: %s", rxn.id, cs)
    logger.debug("%s", compartment)
    reaction_compartment = set()
    for p in cs:
        value = cs[p]
        met_id = p[0]
        if not p[1] in compartment:
            print(rxn.data)
        cmp = compartment[p[1]]
        if cmp is not None and not len(cmp.strip()) == 0:
            met_id += '_' + cmp.strip()
        else:
            cmp = ''
        if met_id not in metabolites:
            metabolites[met_id] = value
            if not p[0] in cpd_cmp:
                cpd_cmp[p[0]] = set()
            cpd_cmp[p[0]].add(cmp)
            reaction_compartment.add(cmp)
        else:
            print('!', rxn.id)

    cmp_token = get_cmp_token(reaction_compartment)
    #print(rxn.id, cmp_token, reaction_compartment)
    crxn['name'] = rxn.id
    if rxn_cmp_suffix is not None and not len(rxn_cmp_suffix.strip()) == 0:
        crxn['id'] = rxn.id + '_' + rxn_cmp_suffix
        crxn['name'] += " [{}]".format(rxn_cmp_suffix)
    else:
        if cmp_token is not None and len(cmp_token) > 0:
            crxn['id'] = rxn.id + '_' + cmp_token
            crxn['name'] += " [{}]".format(cmp_token)
        else:
            crxn['id'] = rxn.id
            
    rev = rxn.data[rev]  # direction
    if lb is None and ub is None:

        lb = -1000
        ub = 1000
        # print('reversibility', rxn.data['reversibility'], 'direction', rxn.data['direction'])
        if rev == '>':
            lb = 0
        elif rev == '<':
            ub = 0
    crxn['lower_bound'] = lb
    crxn['upper_bound'] = ub
    crxn['metabolites'] = metabolites
    crxn['annotation'] = {
        'seed.reaction': rxn.id,
        'seed.compartment': ';'.join(map(lambda x: x[0] + ':' + x[1], compartment.items()))
    }

    logger.debug("%s: %s", rxn.id, metabolites)

    return crxn, cpd_cmp


def modelseed_to_cobra_compound(cpd, cmp='z'):
    ccpd = {
        'id': 'cpd',
        'name': 'Compound',
        'compartment': 'c',
        'charge': 0,
        'formula': '',
        'annotation': {}
    }
    if cmp is not None and not len(cmp.strip()) == 0:
        ccpd['id'] = cpd.id + '_' + cmp
    else:
        ccpd['id'] = cpd.id
    ccpd['name'] = cpd.data['name']
    if cmp is not None and len(cmp) > 0:
        ccpd['name'] += " [{}]".format(cmp)
    ccpd['formula'] = '' if cpd.data['formula'] is None or not type(cpd.data['formula']) == str else cpd.data['formula']
    ccpd['charge'] = cpd.data['charge']
    ccpd['compartment'] = cmp
    ccpd['annotation'] = {'seed.reaction': cpd.id}
    return ccpd


class EscherModelSEEDFactory:
    
    def __init__(self, ms):
        self.modelseed = ms
        self.cmp_config = {}
    
    def build(self):
        all_compounds = []
        all_reactions = []
        all_reaction_compounds = {}
        for rxn_id in self.modelseed.reactions:
            rxn = self.modelseed.get_seed_reaction(rxn_id)
            if not rxn.is_obsolete:
                cmps = set(map(lambda x : x[1], rxn.cstoichiometry.keys()))
                if len(cmps - set(self.cmp_config)) == 0:
                    cobra_reaction, cpd_cmp = modelseed_to_cobra_reaction(rxn, self.cmp_config, None)
                    for cpd_id in cpd_cmp:
                        if cpd_id not in all_reaction_compounds:
                            all_reaction_compounds[cpd_id] = set()
                        for cmp_id in cpd_cmp[cpd_id]:
                            all_reaction_compounds[cpd_id].add(cmp_id)
                    all_reactions.append(cobra_reaction)

        for cpd_id in all_reaction_compounds:
            for cmp_id in all_reaction_compounds[cpd_id]:
                cpd = self.modelseed.get_seed_compound(cpd_id)
                all_compounds.append(modelseed_to_cobra_compound(cpd, cmp_id))

        base_model = {
            'id': 'modelseed',
            'version': '1',
            'compartments': {'z': 'any'},
            'metabolites': all_compounds,
            'reactions': all_reactions,
            'genes': [],
            'annotation': {}
        }
        return base_model
