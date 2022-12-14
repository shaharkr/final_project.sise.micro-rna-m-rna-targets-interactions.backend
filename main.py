from flask import Flask, jsonify, send_from_directory, request
from flask_cors import CORS, cross_origin
import os
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

DB_URL = 'postgresql://postgres:postgrespass@132.73.84.177:8080/micro-message-rna-project'

app.config['SQLALCHEMY_DATABASE_URI'] = DB_URL
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # silence the deprecation warning

db = SQLAlchemy(app)

class Organisim(db.Model):
    __tablename__ = 'organisms'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(200), nullable=False)

    def __init__(self, id, name):
        self.id = id
        self.name = name
    
    def __repr__(self):
        return f"<Organism {self.name}>"

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


@app.route('/api/organisms', methods=['GET'])
def handle_cars():
    cars = Organisim.query.all()
    results = [
        {
            "id": car.id,
            "name": car.name
        } for car in cars]

    return {"count": len(results), "orgs": results}

@app.route('/api/interaction/<int:int_id>', methods=['GET'])
def abc(int_id):
    cars = Interaction.query.filter_by(data_set_id=int_id).all()
    results = [
        {
            "data_set_id": car.data_set_id,
            "index": car.index,
            "sequence": car.sequence
        } for car in cars]

    return {"count": len(results), "inters": results}


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
