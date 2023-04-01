from flask import Flask, jsonify, request
from flask_cors import CORS, cross_origin
from dal.db_connection import init_db_connector
import dal.organisms as organisms
from flask_compress import Compress
import dal.interactions as interactions
from configurator import Configurator
from waitress import serve


app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

Compress(app)
app.config["COMPRESS_REGISTER"] = False
app.config['COMPRESS_ALGORITHM'] = 'gzip'


compress = Compress()
compress.init_app(app)

init_db_connector(app)

@app.route('/api')
@cross_origin("*")
def hello_world():
    to_ret = {
        'status': 'ok',
        'value': 'Hello from micro-message RNA project!@'
    }
    return jsonify(to_ret)

@app.route('/')
@cross_origin("*")
def end_without_api():
    to_ret = {
        'status': 'ok',
        'value': 'Go to api'
    }
    return jsonify(to_ret)

@app.route('/api/organisms/details', methods=['GET'])
@compress.compressed()
def get_organisms_details():
    organisms_list = []
    try:
        with_search_options = request.args.get('searchOptions')
        with_search_options = eval(with_search_options.lower().capitalize())
    except Exception as e:
        print('searchOptions is null or not set properly (as true or flase)')
        with_search_options = False
    try:
        organisms_list = organisms.get_organisms(with_options=with_search_options)
    except Exception as e:
        print(f'app failed to get organisms. error: {str(e)}')
    return organisms_list


@app.route('/api/organisms/datasets/<int:data_set_id>/interactions', methods=['GET'])
def get_data_set_interactions(data_set_id):
    interactions = []
    try:
        interactions = organisms.get_data_set_interactions(data_set_id)
    except Exception as e:
        print(f'app failed to get interactions of data set id- {data_set_id}. error: {str(e)}')
    return interactions


@app.route('/api/interactions', methods=['GET'])
def get_interactions():
    interactions_result = []
    try:
        data_sets_ids = request.args.getlist('datasetsIds')
        seed_families = request.args.getlist('seedFamilies')
        mirna_ids = request.args.getlist('miRnaIds')
        mirna_seqs = request.args.getlist('miRnaSeqs')
        sites = request.args.getlist('sites')
        gene_ids = request.args.getlist('geneIds')
        regions = request.args.getlist('regions')
        interactions_result = interactions.get_interactions(data_sets_ids, seed_families, mirna_ids, mirna_seqs, sites, gene_ids, regions)
    except Exception as e:
        print(f'app failed to get interactions. error: {str(e)}')
    return interactions_result


@app.route('/api/generalSearchInteractions/interactions', methods=['GET'])
def get_general_interactions():
    interactions_result = []
    try:
        query_string = request.args.get('q')
        interactions_result = interactions.get_interactions_gneral_search(query_string)
    except Exception as e:
        print(f'app failed to get general interactions. error: {str(e)}')
    return interactions_result

@app.route('/api/organisms/datasets/<int:data_set_id>/interactions/download', methods=['GET'])
def get_full_data_set(data_set_id):
    response = None
    try:
        response = interactions.download_data(data_set_id, path=None)
    except Exception as e:
        print(f'app failed to get general interactions. error: {str(e)}')
    return response
    

if __name__ == '__main__':
    confg = Configurator()
    mode = confg.get_mode()
    if mode == 'dev':
        app.run(debug=True, host='0.0.0.0', port=5000)
    else:
        serve(app, host='0.0.0.0', port=5000, threads=50)


