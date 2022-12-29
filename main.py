from flask import Flask, jsonify, send_from_directory, request
from flask_cors import CORS, cross_origin
from dal.db_connection import init_db_connector
import dal.organisms as organisms

app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

init_db_connector(app)


@app.route('/api')
@cross_origin("*")
def hello_world():
    to_ret = {
        'status': 'ok',
        'value': 'Hello from micro-message RNA project!'
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
def get_organisms_details():
    organisms_list = []
    try:
        organisms_list = organisms.get_organisms(with_options=True)
    except Exception as e:
        print(f'failed to get organisms. error: {str(e)}')
    return organisms_list


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

