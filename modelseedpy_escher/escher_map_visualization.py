import os
import cobra
from cobra import Model, Reaction, Metabolite
import cobra.test  #needed for the model, not necessary

os.environ["HOME"] = "C:\\Users\\Jason\\Documents\\Agronne"  # needed this for cobrakbase import to work properly
import cobrakbase
from cobrakbase.core.fbahelper import FBAHelper  #fbahelper is what you need from cobrakbase

#jinja / display imports
import jinja2
from jinja2 import Environment, PackageLoader, select_autoescape
from IPython.core.display import HTML
from IPython.display import IFrame
from IPython.core.display import display

#import json to read and write files
import json

# generate the api
from modelseed_escher.core.eschermapapi import EscherMapAPIBiGG as EscherMapAPI  #need escher map API
api = EscherMapAPI()


def hash_model_reactions(model, helper_function=FBAHelper.modelseed_id_from_cobra_reaction):
    # create a hash linking modelseed_id: cobra_reaction_object
    hsh = {}
    for rxnobj in model.reactions:
        if helper_function is None:
            hsh[rxnobj.id] = rxnobj
        else:
            hsh[helper_function(rxnobj)] = rxnobj

    return hsh


def find_reaction_overlap_and_genes(escher_map, model_hash, return_reaction_list=False):
    # HELPER FUNCTION - CALLED IN GENERATE TABLE
    # look at map modelseed reaction ids
    # look through modelseed ids in the hash from the cobra model
    # find the number of overlaps
    # result is a dictionary of all overlaps linking modelseed_id --> cobra reaction object
    # look through the newly created dict of rections
    # compile a dict of genes that are in the reactions (dict so overlaps ar enot counted, could be a set)
    # return the map id, the ratio of overlap_length / model_hash length, and then gene dict length
    genes = {}
    rxns = {}
    for rxn in escher_map['reactions']:
        if rxn in model_hash:
            rxns[rxn] = model_hash[rxn]
    for rxn in rxns:
        genes[rxns[rxn].genes] = rxn

    if return_reaction_list:
        return escher_map['id'], 100 * len(rxns) / len(model_hash), len(genes), rxns
    return escher_map['id'], 100 * len(rxns) / len(model_hash), len(genes)


def generate_table(model_hash, maps, return_reaction_list=False):
    # for each map, create its own dict
    # create final list of dicts to be used in the html input
    # additionally, create a second dictionary for reaction data, paired map_id : dict of reaction data
    table_lst = []
    reaction_data = {}

    if return_reaction_list:

        for E_map in maps:
            emid, rxn_pct, gene, rxn_dict = find_reaction_overlap_and_genes(E_map, model_hash, return_reaction_list)

            reaction_data[emid] = rxn_dict

            row = {}
            row['id'] = emid
            row['reactions'] = round(rxn_pct, 2)
            row['genes'] = gene
            table_lst.append(row)

        return table_lst, reaction_data

    else:
        for E_map in maps:
            emid, rxn_pct, gene = find_reaction_overlap_and_genes(E_map, model_hash, return_reaction_list)

            row = {}
            row['id'] = emid
            row['reactions'] = round(rxn_pct, 2)
            row['genes'] = gene
            table_lst.append(row)

        return table_lst


def render_html(map_lst, template_dir="", html_filename="maps_out.html"):
    # render the html and display in an iframe
    txt = json.dumps(api.get_map('e_coli_core.Core metabolism').escher_data)
    txt = txt.replace('\n', '\\n')

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        autoescape=jinja2.select_autoescape(['html', 'xml']))
    # Return string of html
    template_string = env.get_template('Template.html').render(inp_data=map_lst, map_data_json_b64=txt)
    html_file = open(html_filename, "w")
    html_file.write(template_string)
    html_file.close()
    display(IFrame(html_filename, '100%', '500px'))


def run_application(cobra_model=cobra.test.create_test_model("textbook"),
                    helper_function=FBAHelper.modelseed_id_from_cobra_reaction,
                    map_list=EscherMapAPI().list_maps(),
                    return_reaction_list=False,
                    template_dir="",
                    html_filename="maps_out.html"):
    model_hash = hash_model_reactions(cobra_model, helper_function)
    map_lst = generate_table(model_hash, map_list, return_reaction_list)
    render_html(map_lst, template_dir)
