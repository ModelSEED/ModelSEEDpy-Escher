import logging
import re

logger = logging.getLogger(__name__)


class EschermapError(Exception):
    """Error in escher map API"""
    pass


class EscherMapAPI:
    def __init__(self, api, default_workspace=93991):
        self.api = api
        self.default_workspace = default_workspace
    
    def list_maps(self,workspaces = None):
        if workspaces is None:
            workspaces = [self.default_workspace]
        input = {
            "workspaces":[],"ids":[],'includeMetadata':1,"type":"KBaseFBA.EscherMap"
        }
        for workspace in workspaces:
            if isinstance(workspace, int):
                input["ids"].append(workspace)
            else:
                input["workspaces"].append(workspace)
        maplist = self.api.ws_client.list_objects(input)
        output = []
        for item in maplist:
            name = None
            reactions = None
            type = None
            if "name" in item[10]:
                name = item[10]["name"]
            if "reactions" in item[10]:
                reactions = item[10]["reactions"].split("|")
            if "type" in item[10]:
                type = item[10]["type"]
            if "description" in item[10]:
                description = item[10]["description"]
            newmap = MSEscherMap(item[1],name,description,type,reactions,item[7])
            output.append(newmap)
        return output
    
    def save_map(self,map,id = None,workspace = None):
        if workspace == None:
            workspace = self.default_workspace
        if map.data == None:
            raise EschermapError("Cannot save map without data!")
        map.set_attributes_from_data()
        map.workspace = self.default_workspace
        if id != None:
            map.id = id
            map.data[0]["map_id"] = id
        kbdata = {"metadata":map.data[0],"layout":map.data[1]}
        input = {
            "objects":[{
                "name":map.id,
                "data":kbdata,
                "type":"KBaseFBA.EscherMap",
                "meta":{"type":map.type,"reactions":"|".join(map.reactions),"name":map.name,"description":map.description},
                "provenance":[{
                    'description': 'cobrakbase.core.eschermapapi:',
                    'input_ws_objects': [],
                    'method': 'save_map',
                    'script_command_line': "",
                    'method_params': [{'workspace': workspace,'id': map.id}],
                    'service': 'cobrakbase.core.eschermapapi',
                    'service_ver': "1.0",
                    # 'time': '2015-12-15T22:58:55+0000'
                }]
            }]
        }
        if isinstance(workspace, int):
            input["id"] = workspace
        else:
            input["workspace"] = workspace
        self.api.ws_client.save_objects(input)
        
    def get_maps(self,ids):
        args = {"objects":[]}
        for id in ids:
            if len(id.split("/")) < 2:
                if isinstance(self.default_workspace, int):
                    id = str(self.default_workspace) + "/" + id
                else:
                    id = self.default_workspace + "/" + id
            args["objects"].append({"ref":id})
        output = self.api.get_objects2(args)
        maps = []
        for data in output["data"]:
            newmap = MSEscherMap(data["info"][1])
            corrected_map = [data["data"]["metadata"],data["data"]["layout"]]
            newmap.set_attributes_from_data(corrected_map)
            newmap.workspace = data["info"][7]
            if len(data["info"]) >= 11 and "type" in data["info"][10]:
                newmap.type = data["info"][10]["type"]
            maps.append(newmap)
        return maps

class MSEscherMap:
    def __init__(self,id,name = None,description = None,type = None, reactions = None,workspace = None,data = None):
        self.id = id
        self.name = name
        self.description = description
        self.type = type
        self.reactions = reactions
        self.workspace = workspace
        self.data = data

    def retreive_data(self):
        info = EscherMapAPI.get_maps([self.id],1)
        self.name = info[0].name
        self.type = info[0].type
        self.reactions = info[0].reactions
        self.data = info[0].data
        
    def set_attributes_from_data(self,data = None):
        if data != None:
            self.data = data
        if "map_name" in self.data[0]:
            self.name = self.data[0]["map_name"]
        if "map_description" in self.data[0]:
            self.description = self.data[0]["map_description"]
        self.reactions = []
        for escherid in self.data[1]["reactions"]:
            rxndata = self.data[1]["reactions"][escherid]
            self.reactions.append(rxndata["bigg_id"])
        self.reactions = list(dict.fromkeys(self.reactions))
          