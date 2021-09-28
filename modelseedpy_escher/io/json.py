import json
from modelseedpy_escher.core import EscherMap


def load_json_map(filename):
    """
    """
    with open(filename, "r") as fh:
        return EscherMap(json.load(fh))


def save_json_map(em, filename):
    with open(filename, "w") as fh:
        fh.write(json.dumps(to_json(em)))


def to_json(em):
    """
    Wrapper for later
    """
    return em.escher_data


def from_json(d):
    """
    Load a escher map from a JSON document.

    :param d:
    :return:
    """
    return EscherMap(d)
