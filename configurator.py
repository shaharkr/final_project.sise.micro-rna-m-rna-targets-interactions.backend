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
        return self.paeser['DATASETS']['path_prefix_of_dataset_location']