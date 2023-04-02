from configparser import ConfigParser
import json

class Configurator:
    def __init__(self):
        self.parser = ConfigParser()
        self.parser.read('config.ini')            
    
    def get_db_url(self):
        host = self.parser['DATABASE']['host']
        port = self.parser['DATABASE']['port']
        user = self.parser['DATABASE']['user']
        password = self.parser['DATABASE']['password']
        db_name = self.parser['DATABASE']['database']
        db_url = f'postgresql://{user}:{password}@{host}:{port}/{db_name}'
        return db_url

    def get_mode(self):
        return self.parser['ENV']['mode']
    
    def get_ensambl_url(self):
        return self.parser['PARAMETERS']['ensambl_url']
    
    def get_ensamble_organisms_names_dict(self):
        dict_str = self.parser['PARAMETERS']['ensambl_orgs_names_dict']
        to_ret = json.loads(dict_str)
        return to_ret