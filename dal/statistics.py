from dal.db_connection import GeneralStats, Interaction, DataSet, Organism
from dal.db_connection import db, cache
from sqlalchemy import func, distinct, text
import pandas as pd
from configurator import Configurator


@cache.memoize()
def update_general_stats():
    num_of_organisms = db.session.query(Organism).count()
    num_od_datasets = db.session.query(DataSet).count()
    num_of_mirna = db.session.query(func.count(distinct(Interaction.mirna_id))).scalar()
    gene_ids = db.session.query(Interaction.Gene_ID).all()
    df = pd.DataFrame(gene_ids, columns=['Gene_ID'])
    df['Gene_ID_sub'] = df['Gene_ID'].apply(lambda x: x.split('|')[0])
    num_of_mrna = len(df['Gene_ID_sub'].unique())
    num_of_interactions = db.session.query(Interaction).count()
    num_of_3utr_interactions = db.session.query(Interaction).filter(Interaction.region == "3utr").count()
    general_stats = GeneralStats(num_of_organisms=num_of_organisms,
                              num_of_datasets=num_od_datasets,
                              num_of_mirna=num_of_mirna,
                              num_of_mrna=num_of_mrna,
                              num_of_interactions=num_of_interactions,
                              num_of_3utr_interactions=num_of_3utr_interactions)
    db.session.add(general_stats)
    db.session.commit()
    if db.session.query(GeneralStats).count() > 1:
        row_to_delete = GeneralStats.query.first()
        db.session.delete(row_to_delete)
        db.session.commit()


@cache.memoize()
def get_general_stats():
    try:
        stats = db.session.query(GeneralStats).all()
        num_of_organisms = stats[0].num_of_organisms
        num_of_datasets = stats[0].num_of_datasets
        num_of_mirna = stats[0].num_of_mirna
        num_of_mrna = stats[0].num_of_mrna
        num_of_interactions = stats[0].num_of_interactions
        num_of_3utr_interactions = stats[0].num_of_3utr_interactions
        result = {"featureName": "general_stats",
                "statistics": {
                    "1Organisms": num_of_organisms,
                    "2Datasets": num_of_datasets,
                    "3miRNA": num_of_mirna,
                    "4mRNA": num_of_mrna,
                    "5Interactions": num_of_interactions,
                    "63utr Interactions": num_of_3utr_interactions
                }}
    except Exception as e:
        print(f'dal failed to get stats. error: {str(e)}')
    return result


def get_query_string_for_one_d(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                               site_types, gene_ids, regions, feature_name):
    where_cond_string = get_where_query_part(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                               site_types, gene_ids, regions)
    new_feature_name = f'"{feature_name}"'
    group_by_string = f" GROUP BY {new_feature_name}"
    select_string = f"SELECT {new_feature_name} as feature_values, count({new_feature_name}) FROM mirna_mrna_interactions"
    q = f"{select_string}{where_cond_string}{group_by_string}"
    return q

def get_where_query_part(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                               site_types, gene_ids, regions):
    seed_families = [f"'{x}'" for x in seed_families]
    mirna_ids = [f"'{x}'" for x in mirna_ids]
    mirna_seqs = [f"'{x}'" for x in mirna_seqs]
    gene_ids = [f"'{x}'" for x in gene_ids]
    regions = [f"'{x}'" for x in regions]
    where_cond_list = []
    data_sets_ids_string_cond = f"data_set_id in ({','.join(data_sets_ids)})" if len(data_sets_ids) > 0 else ""
    where_cond_list.append(data_sets_ids_string_cond) 
    seed_families_string_cond = f"seed_family in ({','.join(seed_families)})" if len(seed_families) > 0 else ""
    where_cond_list.append(seed_families_string_cond)
    mirna_ids_string_cond = f"mirna_id in ({','.join(mirna_ids)})" if len(mirna_ids) > 0 else ""
    where_cond_list.append(mirna_ids_string_cond)
    mirna_seqs_string_cond = f"mirna_sequence in ({','.join(mirna_seqs)})" if len(mirna_seqs) > 0 else ""
    where_cond_list.append(mirna_seqs_string_cond)
    char = '"'
    gene_ids_string_cond = f"{char}Gene_ID{char} in ({','.join(gene_ids)})" if len(gene_ids) > 0 else ""
    where_cond_list.append(gene_ids_string_cond)
    regions_string_cond = f"region in ({','.join(regions)})" if len(regions) > 0 else ""
    where_cond_list.append(regions_string_cond)
    site_types_list = []
    if 'canonical' in site_types:
        site_types_list.append('"Seed_match_canonical" IS TRUE AND "Seed_match_noncanonical" IS FALSE')
    if 'noncanonical' in site_types:
        site_types_list.append('"Seed_match_canonical" IS FALSE AND "Seed_match_noncanonical" IS TRUE')
    if 'other' in site_types:
        site_types_list.append('"Seed_match_canonical" IS FALSE AND "Seed_match_noncanonical" IS FALSE')
    site_types_string_cond = ""
    if len(site_types_list) > 0 and len(site_types_list) != 3:
        site_types_string_cond = f"(({') OR ('.join(site_types_list)}))" if len(site_types_list) > 1 else site_types_list[0]
    where_cond_list.append(site_types_string_cond)    
    where_cond_list = list(filter(lambda x: x != "", where_cond_list))
    where_cond_string = f" WHERE {' AND '.join(where_cond_list)}" if len(where_cond_list) > 0 else ""
    return where_cond_string


def get_top_n_dict(statistics_dict, n):
    if len(statistics_dict) < n:
        return statistics_dict
    sorted_dict = dict(sorted(statistics_dict.items(), key=lambda x: x[1], reverse=True))  # Sort the dictionary by values in descending order
    top_dict = dict(list(sorted_dict.items())[:n])  # Get the top 20 elements and their corresponding values
    other_sum = sum(list(sorted_dict.values())[n:])  # Calculate the sum of the values for the remaining elements
    top_dict['others'] = other_sum  # Add the 'remaining_sum' key to the dictionary
    return top_dict


@cache.memoize()
def get_one_d(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
              site_types, gene_ids, regions, feature_name, text_query=None):
    if text_query is None:
        text_query = get_query_string_for_one_d(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                                                site_types, gene_ids, regions, feature_name)
    query = text(text_query)
    try:
        result = db.session.execute(query)
        statistics = {}
        for row in result:
            statistics[row.feature_values] = row.count
        n = sum(list(statistics.values()))
        for k, v in statistics.items():
            statistics[k] = v/n  # convert count freq to %
        
        # get top 20
        top_11_dict = get_top_n_dict(statistics_dict=statistics, n=11)
        
        # convert keys to strings
        to_ret_dict = {}
        for k, v in top_11_dict.items():
            to_ret_dict[str(k)] = v
        data = {"featureName": feature_name, "statistics": to_ret_dict}
        return data
    except Exception as e:
        print(f'dal failed to get oneD statistics for {feature_name}. error: {str(e)}')


@cache.memoize()
def get_dataset_statistics(dataset_id):
    try:
        main_features_name_lst = Configurator().get_main_features_names()
        main_features_back_to_front_names = Configurator().get_main_features_back_to_front_names()
        data =[]
        for feature_name in main_features_name_lst:
            text_q = None
            if feature_name == 'Gene_ID':
                regex = "'([^|]+)'"
                text_q = f'SELECT t."col" as feature_values, count(t."col") FROM (SELECT substring("{feature_name}" FROM {regex}) AS col FROM mirna_mrna_interactions WHERE data_set_id={dataset_id}) t GROUP BY t."col"'
            feature_dist = get_one_d([str(dataset_id)], [], [], [], [], [], [], feature_name, text_q)
            feature_dist["featureName"] = main_features_back_to_front_names[feature_dist["featureName"]]
            data.append(feature_dist)
        return data
    except Exception as e:
        print(f'dal failed to get dataset statistics for id {dataset_id}. error: {str(e)}')


def get_query_string_for_two_d(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                               site_types, gene_ids, regions, feature_name_1, feature_name_2):
    where_cond_string = get_where_query_part(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                                             site_types, gene_ids, regions)
    new_feature_name_1 = f'"{feature_name_1}"'
    new_feature_name_2 = f'"{feature_name_2}"'
    group_by_string = f" GROUP BY {new_feature_name_1}, {new_feature_name_2}"
    select_string = f"SELECT {new_feature_name_1} as feature_1_value, {new_feature_name_2} as feature_2_value, count(*) FROM mirna_mrna_interactions"
    q = f"{select_string}{where_cond_string}{group_by_string}"
    return q


@cache.memoize()
def get_two_d(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
              site_types, gene_ids, regions, feature1, feature2):
    text_query = get_query_string_for_two_d(data_sets_ids, seed_families, mirna_ids, mirna_seqs, 
                                            site_types, gene_ids, regions, feature1, feature2)
    query = text(text_query)
    try:
        result = db.session.execute(query)
        x = []
        y = []
        for row in result:
            x.append(row.feature_1_value)
            y.append(row.feature_2_value)
        
        data = {"featureX": feature1,
                "secondFeature": feature2,
                "statistics": {
                    "x": x,
                    "y": y
                    }
                }
        return data
    except Exception as e:
        print(f'dal failed to get twoD statistics for {feature1}~{feature2}. error: {str(e)}')

