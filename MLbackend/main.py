import pandas as pd
from full_pipline_positive import full_pipline as full_pipline_positive
from full_pipline_negative import full_pipline as full_pipline_negative
import pickle
from consts.global_consts import list_of_features_for_web as features
from flask import Flask, request
import subprocess
import xgboost as xgb



app = Flask(__name__)


@app.route('/wsl_api')
def hi():
	return 'hello from mrna-mirna wsl'


@app.route('/wsl_api/prediction')
def prediction():
    try:
        miRNA_seq = request.args.get('miRNA_seq')
        mRNA_seq = request.args.get('mRNA_seq')
        site_seq = request.args.get('site_seq')
        org_name = request.args.get('org_name')
        prob = predict_sequences(miRNA_seq, mRNA_seq, site_seq, org_name)
        prob_float = float(prob[0])
        response_data = {'prob': prob_float}
    except Exception as e:
        print("Error with prediction: " + str(e))
    return response_data

if __name__ == '__main__':
    command = [
        "gunicorn",
        "-w", "8",
        "-b", "0.0.0.0:3030",
        "main:app"
    ]
    subprocess.run(command)

def load_model(org_name):
    model = xgb.Booster()
    model.load_model(org_name)
    return model

def extract_features(miRNA=None, mRNA=None, site=None):
    columns = ['key', 'paper name', 'organism', 'miRNA ID', 'miRNA sequence', 'site', 'region', 'valid_row',
                        'full_mrna', 'Gene_ID', 'region count']
    values = [None] * len(columns)
    values[4] = miRNA
    values[5] = site
    values[8] = mRNA
    dataset = pd.DataFrame([values], columns = columns)

    if site != None:
        result_df = full_pipline_positive(dataset)
    else:
        result_df = full_pipline_negative(dataset)
        result_df.drop('site_x', axis=1,inplace=True)

    return result_df


def extract_features_from_sequences(df):

    if df['site'].iloc[0] != None:

        result_df = full_pipline_positive(df)
    else:
        result_df = full_pipline_negative(df)
        result_df.drop('site_x', axis=1,inplace=True)

    return result_df


def get_prediction(seq_features,org_name):

    # load the model (the names ends with .model)
    model = load_model(org_name)
    # leave only the columns that the model trained on (list in consts.global_consts)
    columns_to_keep = list(set(seq_features.columns) & set(features))
    # get the features (requires using the viennaRNA module)
    seq_features = seq_features[columns_to_keep]
    # predict 1 for True 0 for False
    dmatrix = xgb.DMatrix(seq_features)
    pred = model.predict(dmatrix)
    return pred


# def filter_columns(seq_features):
#
#     features_to_remove = ['miRNA ID', 'miRNA sequence', 'Gene_ID', 'new_key', 'MiRBase ID',
#                           'Node of origin (locus)', 'Node of origin (family)', 'Family Age_Group',
#                           'Family Age', 'Locus Age', 'seed_family', 'site', 'region', 'paper region',
#                           'sequence', 'Gene_name', 'mrna_bulge', 'mrna_inter', 'mir_inter',
#                           'mir_bulge', 'duplex_method'
#                           ]
#
#     # removing features
#     features_to_remove = list(set(seq_features.columns).intersection(features_to_remove))
#     seq_features = seq_features.drop(columns=features_to_remove)
#
#     # Remove datatypes that are not numeric
#     seq_features = seq_features.select_dtypes(include=['float64', 'int64'])
#
#     seq_features.drop(['index', 'key'], axis=1)
#     return seq_features


def predict_sequences(miRNA_seq=None, mRNA_seq=None, site_seq=None,org_name="Human.model"):
    organisem_name = {"Human" : "Human.pkl",
                      "Mouse" : "Mouse.model",
                      "Cattle" : "Cattle.model",
                      "Worm" : "Worm.model",
                      "Fly" : "Fly.model"}

    seq_features = extract_features(miRNA_seq, mRNA_seq, site_seq)
    prediction = get_prediction(seq_features, organisem_name[org_name])
    return list(prediction)


# if __name__ == '__main__':
#     print(predict_sequences('UGGAGUGUGACAAUGGUGUUUG',
#                             'CAAUAGAUGUGAGUUAAACUUUAGGAAAAAGGAUUCCCUUUUUUUAAAAAAAAUCAAUACCUCAAAAGCAGGCUUUGGGACAAGAAAACCCCAAAGUGGCCUGCUUUUCCCAUCCCAGGAGCUCAUUAUCCAGUCUGUGCCAACUGAAGUAGGAGACUGACUGUGAGUGCUGGCUAAAAGCCCUGGGUGGUGAGGCUCACAGUACUGGUUUCCAGGAGGAAGAGCCUUUGUGCAUUUGACUGAGGCCAGUUUCUAUGAAGAGCAAGUAGCUGAGGAGAGGUCGAAUUUACUGCUUUUUCCAGGACAAUUCUGGAAGUAAAGAAAAUGUAAUUCAAGCUGGUUAGCUUAAUUUUGUGCCAUUCUUUAACAUAAGAGUAAGCUCUAUUAUGAAAUACAACUUUAAAAAAUUUUAGCUAUAAAUUAUAUAAAUGAUUUUAAAUUGCUGAGGUUUCCUUAGGCAGCUUAUUUAUUUGUUUACAGUUAGACUAUCUGAGUAAAUGGUUCUUUGUGGACCUAGGCAGUUCCUGACUGUUCCACAUGUAGUACAUUGUACCAAAGUUCUUAAUAAGAAUAUUCCCCACAAUCCUGUUCUCUAAAUGUCAAAUAAAGAUUAUUUUCACUAGAUUCAACUUUACAAAAUUUGUUUUAUAUCUGUUAGAAAAUGUACAGACAUAAGUAUUUUCAGUUGACAAAGCAUCAAACCCAGUUCUGCCUAGUGAUAAGUUUCACCCUAGAGUAUGUAUGUAACGUUUUAGCUUAUCCAUCCUUUCUUGGAGCGCCUCCAUUUCCAUUGAAAGCCAGGCUGGAGCAGGACCCUUUUGGAGUAGUGACUCAGUUGCUUCCAAAGCCCCUGCUAUUGUAUGCAGCGCUGACCUGUACUCUUCUUCCCAGGGGAACUCCUGACGAGCUCUUUUUGCAUAAGGCUGGAAAAAAAACAUAAGUAAUAUCACAAUAUCCAUUCUAAAUAUAAAGAACCUUCCUUUUGGACUGGAGUAAAGCUUACAUGCAAAUUUUAUUCUAGUCAUUGGAUCACAAGGGUAGGAGGAUGCACCCCAAAACCCCUACACAGUCAUCUAGAAAAAUAUGUAAAGGCAUUUUGGUUUAUCAUAGCAAUUCAGAGUGCUACUACCAGUGUCUUAGUUUGUAUGUGGUAUACAACAAGUAUCCUGUCCCAAAGGGCUCCCAAUGAGAAGUGCUGCAUAGUCCAAGCUUACAUGUCUUAUAAACAAGUUCAUAAAUGUAUUUUCUUUUUAUGAGAGUUUGACUAAAACUUAUCAGAAUGUUGUUCUUCAUGAAUUACUACUAUACUAAUUACUAUACUAAUAGUGCUCAAAACAAUAUUUUGAAUAUCCUUAUUGGUGUCAAAUUCUGCCUUUUAAUAAGUAGAUGUGAUCUUCAGUUACUGCCAAAAAUUAUUAGGAGACUCAUUUGAUUAAUAAGGCAAGGAAUCAAACUAAACAUUUAGGAGUAAGUUUCUUUCAUUUUCUUCUGUGGUUCAGUAAAGACUGCAUUUAUAGCAUCACUGGUACAAUAUGUAACUUCCCUUAAAGGUUACUACCAAUAAUUCAAACAUACUGAAAGAAUAUAUUUGAUAUGGUGUAGUCCCACUUCUUAAUUUUAAAAGCAACUACCAUAAAACAGAAUUUUACAUGUCUAGAUCUAUUUGAUUUGAAAUUCAGCAUAAGGCUGGAAACCACACUGGUUUGUUUCGUCAGUAAGUAAAAAGGGCAGAAUUUGCCUUGUUAAAGUUUGGCCCCUAUUGAAAUCAGCCCAUACCUGUAAAGAUGACCUCUUUGCUUCUUCUACAGUCACAUUAGCAAAGGGUUCCCAGAAAAUACCUUUUUCCUGUUUCACACGUUCCACUUUGGCAGCUUCAGUUUCAUCUACAAACCCAGUCUGCCAGGGACCAUGAAAAACCAAGCAAAUAGCAACAUGUUAGCACUCUACUAGAUAUGAAAUGGCCACAUAAUUUAAGUGCUGAGUGUUCAACCAUACUAGGCAAAUUCUGGAAGUCACAUGGACUGAUCUAUAAAUACUCUUAGUAUAAUCUGGAUUAAAUCACUCUACUGUGUUCUUCCUUACAAAAUAGUUAGUUAUAGAGUUGUCUUUGCAGGAAAAAAUGAUAGUCUGCUAACCUUUACUGUAUAAGCUAAGAAACUGGCAACAGCAGUGUAUACACAUGUACUUAAAAUCCAGCUUAGAUGUAAUAUAAGUAGGCCAAGUGUGGUGGCUCAUGCCUCUAAUCCCAGGACUUUGGUGGGAGGAAUGCUUGAGCCCAGGAGUUUGAGGUUAUAGCGAGCUAUGAUCACACCACCACACUCCAAUCUGGGCGACAGAGCGAGACUCUGUCUAAAAAAUAAAGUUCAGUGAGUCAAGUGGGCCCAUUUCACCUAAAAAAAUAUAUAUGUAAUAUAUGUAAAAUAUAUAACGUGCAUCAUUACAUGUAAUAUAUCAAAUAAUUCAGUUUCUCCUGGUAGUCUAAAUUUGUGGUUUAACAGUUCCACUAGUAUUUCAAUUAACUAUUAUUUCCUCACACUAAGUUCCCAAGUAUGCAUAGACAUAGGAACCUUGUGAUUCAAAAUUUUGGUUUUAAUUGUAAAAACAGGUCUUGGCUGCAAAGAGAAUAAAAAAGCCAUGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGCCCUGAGCUUACCCUGGCUUGGGAGGCUGCCCAAAUAAAAAAAAAAAAAAAAAAAAAAAAAAGCCAUGCCUAAAUUAUACCCAAAAUAUGGUAGUAAACAUCCUCUCAGGCCUACCUAAUUAUGUGAUCAAGUGUAUUUAAUUAUGAUUAUUAAUAGCUUAGGGAUCCUCAUUAGUCAUCUCACUGCCUGGAUAAUCAGUGCUAUUACACCCCAAACUACAGGCAGAGGAACUACAAGUGUCCUUUUUUGAGACAGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGGCAUGAUCUUGGCUUACUGCAAUCUCUACCUCCCCAGUUCAAGCAAUUCUCCUGCCUCAGCCUCCCGAGUAGCUGGGACUGCAGGCGCGUACCACCAUGCCCGGCUAAUUUUUGUAUUUUUAGUAGAAACAGGGUUUCCCCAUGUUGGCCAGGCUGGUCUUGAACUCCUGACCUCAUGAUCCACCCACCUUGGCCUCCCAAAGUGCUGGGAUUACAGGCAUGAGCCGCCGCGCCUGGCCCCAAGUGUUUUACUUCUGUGGCCUCAAUUGCUUAGUAAAAGUCGUCAGCACUGUGCUCAGCACUGUGAAAAGUUUAAGCUAAGAAGACACUCAGGUUGGUUACAGUUACAUAAUUUCGGAAAAAAAAAAAUAGGACAAAGAAUAAAACAUGAGCCAUGUAAAGAAAGUACCUUCUCCAGAAAGGAUACAACUUUGUUCUUAGUUUCUUCAGAACUGAGACACUGUGGAACUAUCUCUUCAUGAUUUGUUCCCUAUUAAAAAAUUGAUGAACAAGACCAAUUUUAACAUUUCAACAACUUGCAAAGACAAAUAGAUGAGUUUUAUGAGCUACCUUAACCAAAAAUUCAUUUUAAGUAAGAGUCCCAGAGGAUCCUCAAAGGUGAUAAACUCAUGAUUCCUUCAGGGUCCCUAAGGAUAAUACAAAAUUAACUUCUGUACAAGUGUUGUAAGCUUUAAUUACUUCUGCUGGGUCAUACUAAUGCUUCUGGACUCCCUUAUGAUCAUAGGCAUGACUAUACAGCAAAUGCAGUAAGAGUAAAGGAGCAACCCAUCUUACAGGUUUAGGUUCAUCAGCCUAACCCUUAUGACUGAUAGCACAAAAUGAAAUGUAUUAUCAUUUGACCCAAAAAUACUAUCUGCUGGAAGACUGUGUCUGUGUCUGUCUGUAUCAGUAGGCCUGCUGUGUAUGCCUGUUCUGGUCCUUGUUUAUUCAGAGGACCUUACGAAAUUCACUUCAUUUAUCUAAGCCUCAUUUUGAGAAGCUGUAAAAGAGAUAACGAGUAAUGUACCCUUCAGACAAUUUUCCGAUUGCAAUACAGAAGCAGUUCAAUAAAUGUUUUGGGAUUGUUCUGGAAUAUUUGAAAUAUUAAAAUGGUUUGAAAGUCA',
#                             'CAAAUAGCAACAUGUUAGCACUCUA',
#                             org_name='Human.model'))



