class EscherReaction:

    def __init__(self, reaction_id, name='', reversibility=1):
        self.id = reaction_id
        self.segments = {}
        self.annotation = {}
        self._metabolites = []
        self.genes = []
        self.gene_reaction_rule = ''

    @property
    def bigg_id(self):
        return self.id
