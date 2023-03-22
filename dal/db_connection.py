from flask_sqlalchemy import SQLAlchemy
from configurator import Configurator
from flask_executor import Executor
from flask_caching import Cache

db = SQLAlchemy()
executor = Executor()
cache = Cache(config={'CACHE_TYPE': 'SimpleCache'})

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
        cache.init_app(app)
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
    data_set_mb = db.Column(db.Float, nullable=False)
    interactions_amount = db.Column(db.Integer, nullable=False)
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


class GeneIdOption(db.Model):
    __table__ = db.Table('distincted_gene_id', db.metadata,
        db.Column("Gene_ID", db.String(200), nullable=False, primary_key=True),
        db.Column("data_set_id", db.Integer, primary_key=True)
        )
    
    def __repr__(self):
        return f"<GeneIdOption: Gene_ID-{self.Gene_ID}, data_set_id-{self.data_set_id}>"


class SeedFamilyOption(db.Model):
    __table__ = db.Table('distincted_seed_family', db.metadata,
        db.Column("seed_family", db.String(200), nullable=False, primary_key=True),
        db.Column("data_set_id", db.Integer, primary_key=True)
        )
    
    def __repr__(self):
        return f"<SeedFamilyOption: seed_family-{self.seed_family}, data_set_id-{self.data_set_id}>"


class mirnaIdOption(db.Model):
    __table__ = db.Table('distincted_mirna_id', db.metadata,
        db.Column("mirna_id", db.String(200), nullable=False, primary_key=True),
        db.Column("data_set_id", db.Integer, primary_key=True)
        )
    
    def __repr__(self):
        return f"<mirnaIdOption: mirna_id-{self.mirna_id}, data_set_id-{self.data_set_id}>"


class Interaction(db.Model):
    __tablename__ = 'mirna_mrna_interactions'
    index = db.Column(db.Integer, primary_key=True)
    data_set_id = db.Column(db.Integer, primary_key=True)
    mirna_id = db.Column(db.String(200))
    mirna_sequence = db.Column(db.String(200))
    seed_family = db.Column(db.String(200))
    site = db.Column(db.String(200))
    region = db.Column(db.String(200))
    start = db.Column(db.Integer)
    end = db.Column(db.Integer)
    mrna_bulge = db.Column(db.String(200))
    mrna_inter = db.Column(db.String(200))
    mir_inter = db.Column(db.String(200))
    mir_bulge = db.Column(db.String(200))
    Gene_ID = db.Column(db.String(200))
    Energy_MEF_Duplex = db.Column(db.Float)
    Seed_match_start = db.Column(db.Float)
    MRNA_Dist_to_start = db.Column(db.Float)
    MRNA_Dist_to_end = db.Column(db.Float)
    organism = db.Column(db.String(200))
    paper_name = db.Column(db.String(200))
    Seed_match_noncanonical = db.Column(db.String(200))
    Seed_match_canonical = db.Column(db.Boolean)
    
    def __repr__(self):
        return f"<Interaction: index-{self.index}, data_set_id-{self.data_set_id}>"
    
