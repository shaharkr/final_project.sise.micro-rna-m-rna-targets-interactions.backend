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