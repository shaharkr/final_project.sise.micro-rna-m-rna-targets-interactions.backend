from flask import Flask, jsonify, send_from_directory
from flask_cors import CORS, cross_origin
import os

app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'


@app.route('/')
@cross_origin("*")
def hello_world():
    return 'hi from micro message rna project'

@app.route('/app')
@cross_origin('*')
def get_app():
    return send_from_directory('./../sise.micro-rna-m-rna-targets-interactions.frontend/react-app/build', 'index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
