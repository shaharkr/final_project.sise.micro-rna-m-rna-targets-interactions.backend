from flask_sqlalchemy import SQLAlchemy
from configurator import Configurator
from flask_executor import Executor

db = SQLAlchemy()
executor = Executor()

def init_db_connector(app):
    DB_URL = None
    try:
        config = Configurator()
                
        DB_URL = config.get_db_url()
    except Exception as e:
        print(f'failed to get the connection string. error: {str(e)}')
    try:
        app.config['SQLALCHEMY_DATABASE_URI'] = DB_URL
        app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # silence the deprecation warning

        # db = SQLAlchemy(app)
        db.init_app(app)
        executor.init_app(app)
    except Exception as e:
        print(f'failed to init the db connection. error: {str(e)}')

class Organism(db.Model):
    __tablename__ = 'organisms'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(200), nullable=False)
    data_sets = db.relationship('DataSet', lazy=True,backref='organism')

    def __init__(self, id, name):
        self.id = id
        self.name = name
    
    def __repr__(self):
        return f"<Organism: id-{self.id}, name-{self.name}>"

class DataSet(db.Model):
    __tablename__ = 'data_set'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(200), nullable=False)
    data_set_mb = db.Column(db.Float)
    interactions_amount = db.Column(db.Integer)
    organism_id = db.Column(db.Integer, db.ForeignKey('organisms.id'))
    
    def __repr__(self):
        return f"<DataSet: id-{self.id}, oranism-{self.organism}>"


class SiteOption(db.Model):
    __table__ = db.Table('distincted_sites', db.metadata,
        db.Column("site", db.String(200), nullable=False, primary_key=True),
        db.Column("data_set_id",db.Integer, primary_key=True)
        )
    def __repr__(self):
        return f"<SiteOption: site-{self.site}, data_set_id-{self.data_set_id}>"


class RegionOption(db.Model):
    __table__ = db.Table('distincted_regions', db.metadata,
        db.Column("region", db.String(200), nullable=False, primary_key=True),
        db.Column("data_set_id", db.Integer, primary_key=True)
        )
    
    def __repr__(self):
        return f"<RegionOption: region-{self.region}, data_set_id-{self.data_set_id}>"
