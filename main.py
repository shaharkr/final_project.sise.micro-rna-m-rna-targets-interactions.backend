from flask import Flask, jsonify, send_from_directory, request
from flask_cors import CORS, cross_origin
import os
from flask_sqlalchemy import SQLAlchemy
from configurator import Configurator


app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

config = Configurator()
        
DB_URL = config.get_db_url()

app.config['SQLALCHEMY_DATABASE_URI'] = DB_URL
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # silence the deprecation warning

db = SQLAlchemy(app)

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

class Interaction(db.Model):
    __tablename__ = 'mirna_mrna_interactions'
    data_set_id = db.Column(db.Integer, primary_key=True)
    index = db.Column(db.Integer, primary_key=True)
    sequence = db.Column(db.String(33000), nullable=True)

    def __init__(self, data_set_id, index, sequence):
        self.data_set_id = data_set_id
        self.index = index
        self.sequence = sequence
    
    def __repr__(self):
        return f"<Interaction {self.name}>"


@app.route('/api')
@cross_origin("*")
def hello_world():
    to_ret = {
        'status': 'ok',
        'value': 'Hello from micro-message RNA project!'
    }
    return jsonify(to_ret)


@app.route('/api/organisms/details', methods=['GET'])
def get_organisms_details():
    data_sets_dict = get_data_sets()
    get_search_options()
    organisms_dict = {}
    for data_set in data_sets_dict.values():
        if data_set["organism"].id not in organisms_dict:
            organisms_dict[data_set["organism"].id] = {"id": data_set["organism"].id,
                                                    "name": data_set["organism"].name,
                                                    "datasets": []}
        organisms_dict[data_set["organism"].id]["datasets"].append(data_set)
        data_set.pop("organism")
    organisms_list = list(organisms_dict.values())
    return organisms_list

def get_data_sets():
    data_sets = db.session.query(DataSet)
    data_sets_dict = {}
    for data_set in data_sets:
        data_sets_dict[data_set.id] = { "id": data_set.id,
                                        "name": data_set.name,
                                        "interactionsAmount": data_set.interactions_amount,
                                        "datasetMB": data_set.data_set_mb,
                                        "searchOptions": {"seedFamilies": [], "miRnaIds": [], "miRnaSeqs": [],
                                                        "sites": [], "geneIds": [], "regions": []},
                                        "organism" : data_set.organism}
    return data_sets_dict

def get_search_options():
    pass

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
