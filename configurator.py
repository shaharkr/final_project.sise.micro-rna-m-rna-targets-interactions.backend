from configparser import ConfigParser
import json

class Configurator:
    def __init__(self):
        self.parser = ConfigParser()
        self.parser.read('config.ini')
        self.statistic_parser = ConfigParser()
        self.statistic_parser.read('statistics_config.ini')
    
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
      
    def get_path_prefix_of_dataset_location(self):
        return self.parser['DATASETS']['path_prefix_of_dataset_location']
    
    def get_path_prefix_to_save_new_csv(self):
        return self.parser['DATASETS']['path_prefix_to_save_new_csv']
    
    def convert_string_to_list_of_dicts(self, string):
        # Convert the string to a list of dictionaries
        try:
            # Parse the string as JSON
            data = json.loads(string)
            
            # Ensure that the parsed data is a list
            if isinstance(data, list):
                # Check if each item in the list is a dictionary
                if all(isinstance(item, dict) for item in data):
                    return data
                else:
                    raise ValueError("The input string does not represent a list of dictionaries.")
            else:
                raise ValueError("The input string does not represent a list.")
        except json.JSONDecodeError:
            raise ValueError("Invalid string format or invalid JSON.")
    
    def get_statistic_features_details(self):
        data = self.statistic_parser['DETAILS']['features_details']
        return self.convert_string_to_list_of_dicts(data)
    
    def get_features_types(self):
        data = self.statistic_parser['DETAILS']['features_types']
        return self.convert_string_to_list_of_dicts(data)
    
    def get_main_features_names(self):
        return eval(self.statistic_parser['FEATURES']['main_features_names_list'])
    
    def get_main_features_back_to_front_names(self):
        dict_str = self.statistic_parser['FEATURES']['main_features_back_to_front_names']
        return json.loads(dict_str)
