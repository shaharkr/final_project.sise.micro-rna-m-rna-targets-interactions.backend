from dal.db_connection import DataSet, db


def get_organisms():
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