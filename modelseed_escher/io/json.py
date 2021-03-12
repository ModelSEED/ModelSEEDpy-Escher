import json
from modelseed_escher.core import EscherMap


def load_json_map(filename):
    """
    """
    with open(filename, "r") as fh:
        return EscherMap(json.loads(fh.read()))
