from dal.db_connection import DataSet, db, executor
from dal.db_connection import SeedFamilyOption, SiteOption, GeneIdOption, RegionOption

def get_organisms(with_options=False):
    data_sets_dict = get_data_sets(with_options)
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

def get_data_sets(with_options=False):
    if with_options:
        search_options_dict = get_search_options()
    data_sets = db.session.query(DataSet)
    data_sets_dict = {}
    for data_set in data_sets:
        data_sets_dict[data_set.id] = { "id": data_set.id,
                                        "name": data_set.name,
                                        "interactionsAmount": data_set.interactions_amount,
                                        "datasetMB": data_set.data_set_mb,
                                        "organism" : data_set.organism}
        if with_options:
            data_sets_dict[data_set.id]["searchOptions"] = {"seedFamilies": None, "miRnaIds": None, "miRnaSeqs": None,
                                                            "sites": None, "geneIds": None, "regions": None}
    if with_options:
        for op_name, options_of_data_sets_dict in search_options_dict.items():
            for data_set_id, options_list in options_of_data_sets_dict.items():
                if data_set_id in data_sets_dict:
                    data_sets_dict[data_set_id]["searchOptions"][op_name] = options_list
    return data_sets_dict

def get_search_options():
    options_details_dict = {"seedFamilies":[SeedFamilyOption, "seed_family"],
                            "sites": [SiteOption, "site"],"geneIds":[GeneIdOption, "Gene_ID"],  "regions": [RegionOption, "region"]}
    workers = []
    for op_name, op_details in options_details_dict.items():
        workers.append(executor.submit(get_search_option, op_name, op_details[0], op_details[1]))
    return dict([w.result() for w in workers])


def get_search_option(name, option, column_name):
    res = {}
    print(f'---start extraction of {name}')
    try:
        data = db.session.query(option).all()
        for d in data:
            if d.data_set_id not in res:
                res[d.data_set_id] = []
            res[d.data_set_id].append(getattr(d, column_name))
        print(f'---end extraction of {name}')
    except Exception as e:
        print(f'---extraction failed- {name}. error: {str(e)}')
    return (name, res)