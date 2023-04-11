from configparser import ConfigParser

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
    
    def get_path_prefix_of_dataset_location(self):
        return self.parser['DATASETS']['path_prefix_of_dataset_location']
    
    def get_path_prefix_to_save_new_csv(self):
        return self.parser['DATASETS']['path_prefix_to_save_new_csv']
    
    # def update_max_value_in_temp_files_folder(self, new num):
    #     if reset_mode is not None:
    #         set_new_number = reset_mode
    #         self.parser.set('DATASETS', 'next_num_for_temp_files_folder', str(set_new_number))
    #         with open('config.ini', 'w') as configfile:
    #             self.parser.write(configfile)
    #     else:
    #         set_new_number = int(self.parser['DATASETS']['next_num_for_temp_files_folder']) + 1
    #         self.parser.set('DATASETS', 'next_num_for_temp_files_folder', str(set_new_number))
    #         with open('config.ini', 'w') as configfile:
    #             self.parser.write(configfile)
    #     return set_new_number
