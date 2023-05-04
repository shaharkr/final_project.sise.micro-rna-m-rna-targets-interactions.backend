# server to run inside a wsl environement

from flask import Flask, request
import RNA

app = Flask(__name__)


@app.route('/wsl_api')
def hi():
        return 'hello from mrna-mirna wsl'


@app.route('/wsl_api/duplexfold')
def duplexfold():
        string1 = request.args.get('string1')
        string2 = request.args.get('string2')
        duplex = RNA.duplexfold(string1,string2)
        return {'structure': str(duplex.structure), 'energy': str(duplex.energy), 'i': str(duplex.i), 'j': str(duplex.j>
if __name__ == '__main__':
    app.run(port=3030)
