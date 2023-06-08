from consts.global_consts import DATA_PATH, ROOT_PATH, DATA_PATH_INTERACTIONS
from utilsfile import read_csv, to_csv
import pandas as pd


def comper():
    dataset_old = DATA_PATH_INTERACTIONS / "comper/old.csv"
    dataset_new = DATA_PATH_INTERACTIONS / "comper/new.csv"
    dataset_comper = DATA_PATH_INTERACTIONS / "comper/comper.csv"
    old = read_csv(dataset_old)
    new = read_csv(dataset_new)

    #df = set(new.columns.values) - set(old.columns.values)
    df = set(old.columns.values) - set(new.columns.values)

    s = list(df)
    df = pd.DataFrame(s)
    #df.to_csv(dataset_comper)
    print(df, "CHANGE")

def comper_intersting():
    dataset_old = DATA_PATH_INTERACTIONS / "comper/old.csv"
    dataset_new = "/sise/home/efrco/efrco-master/data/train/stratify/mockMirna_unambiguous_human_ViennaDuplex_negative_features_train_stratify_method.csv"
    dataset_comper = DATA_PATH_INTERACTIONS / "comper/comper1.csv"
    old = read_csv(dataset_old)
    new = read_csv(dataset_new)

    col_list_new = list(new.columns)
    all_features_new = col_list_new[col_list_new.index("Seed_match_interactions_all"):]
    print(len(all_features_new))

    col_list_old = list(old.columns)
    all_features_old = col_list_old[col_list_old.index("Seed_match_compact_A"):]
    print(len(all_features_old))


    # df = set(all_features_new) - set(all_features_old)
    df = set(all_features_new) - set(all_features_old)


    print(df)
    s = list(df)
    df = pd.DataFrame(s)
    df.to_csv(dataset_comper)

def split_features_to_groups_old():
    dataset_old = DATA_PATH_INTERACTIONS / "comper/old.csv"
    dataset_new = DATA_PATH_INTERACTIONS / "comper/new.csv"
    data = read_csv(dataset_old)
    new = read_csv(dataset_new)
    features = data.columns
    Seed_match_interactions_columns = [x for x in features if 'Seed_match_interactions' in x]
    miRNA_pairing_columns = [x for x in features if 'miRNAMatchPosition' in x]
    miRNA_pairing_count_columns = [x for x in features if 'miRNAPairingCount' in x]
    Hot_pairing_mRNA_columns = [x for x in features if 'HotPairingMRNA' in x]
    Hot_pairing_miRNA_columns = [x for x in features if 'HotPairingMirna' in x]
    Composition_columns = [x for x in features if '_comp' in x and 'compact' not in x]
    Energy_columns = [x for x in features if 'MEF' in x]
    Accessibility_columns = [x for x in features if 'Acc_P' in x]
    Distance_columns = [x for x in features if 'Dist' in x]
    temp_list = Seed_match_interactions_columns + miRNA_pairing_columns + Composition_columns + Energy_columns + Accessibility_columns + Distance_columns + miRNA_pairing_count_columns + Hot_pairing_mRNA_columns + Hot_pairing_miRNA_columns
    general_features_columns = [x for x in features if x not in temp_list]
    Seed_match_columns = [x for x in general_features_columns if 'Seed_match' in x]
    general_features_columns = [x for x in general_features_columns if x not in Seed_match_columns]

    Seed_match = data[Seed_match_columns]
    Seed_match_interactions = data[Seed_match_interactions_columns]
    miRNA_pairing = data[miRNA_pairing_columns]
    miRNA_pairing_count = data[miRNA_pairing_count_columns]
    Hot_pairing_mRNA = data[Hot_pairing_mRNA_columns]
    Hot_pairing_miRNA = data[Hot_pairing_miRNA_columns]
    Composition = data[Composition_columns]
    Energy = data[Energy_columns]
    Accessibility = data[Accessibility_columns]
    Distance = data[Distance_columns]
    general_features = data[general_features_columns]
    # return Seed_match, Seed_match_interactions, miRNA_pairing, miRNA_pairing_count, Hot_pairing_mRNA, \
    #        Hot_pairing_miRNA, Composition, Energy, Accessibility, Distance, general_features
    return general_features.columns


def split_features_to_groups_new():
    dataset_new = DATA_PATH_INTERACTIONS / "comper/new.csv"
    data = read_csv(dataset_new)
    features = data.columns
    Seed_match_interactions_columns = [x for x in features if 'Seed_match_interactions' in x]
    miRNA_pairing_columns = [x for x in features if 'miRNAMatchPosition' in x]
    miRNA_pairing_count_columns = [x for x in features if 'miRNAPairingCount' in x]
    Hot_pairing_mRNA_columns = [x for x in features if 'HotPairingMRNA' in x]
    Hot_pairing_miRNA_columns = [x for x in features if 'HotPairingMirna' in x]
    Composition_columns = [x for x in features if '_comp' in x and 'compact' not in x]
    Energy_columns = [x for x in features if 'MEF' in x]
    Accessibility_columns = [x for x in features if 'Acc_P' in x]
    Distance_columns = [x for x in features if 'Dist' in x]
    temp_list = Seed_match_interactions_columns + miRNA_pairing_columns + Composition_columns + Energy_columns + Accessibility_columns + Distance_columns + miRNA_pairing_count_columns + Hot_pairing_mRNA_columns + Hot_pairing_miRNA_columns
    general_features_columns = [x for x in features if x not in temp_list]
    Seed_match_columns = [x for x in general_features_columns if 'Seed_match' in x]
    general_features_columns = [x for x in general_features_columns if x not in Seed_match_columns]

    Seed_match = data[Seed_match_columns]
    Seed_match_interactions = data[Seed_match_interactions_columns]
    miRNA_pairing = data[miRNA_pairing_columns]
    miRNA_pairing_count = data[miRNA_pairing_count_columns]
    Hot_pairing_mRNA = data[Hot_pairing_mRNA_columns]
    Hot_pairing_miRNA = data[Hot_pairing_miRNA_columns]
    Composition = data[Composition_columns]
    Energy = data[Energy_columns]
    Accessibility = data[Accessibility_columns]
    Distance = data[Distance_columns]
    general_features = data[general_features_columns]
    # return Seed_match, Seed_match_interactions, miRNA_pairing, miRNA_pairing_count, Hot_pairing_mRNA, \
    #        Hot_pairing_miRNA, Composition, Energy, Accessibility, Distance, general_features
    return general_features.columns


def comper_pos_neg():
    dataset_neg ="/sise/home/efrco/efrco-master/data/train/stratify/mockMirna_unambiguous_human_ViennaDuplex_negative_features_train_stratify_method.csv"
    dataset_pos = "/sise/home/efrco/efrco-master/data/positive_interactions/darnell_human_ViennaDuplex_features.csv"
    Q_CLASH="/sise/home/efrco/efrco-master/data/positive_interactions/qclash_melanoma_human_ViennaDuplex_features.csv"
    dataset_comper = DATA_PATH_INTERACTIONS / "comper/comper1.csv"
    old = read_csv(dataset_pos)
    new = read_csv(dataset_neg)
    open = read_csv(Q_CLASH)

    col_list_new = list(new.columns)
    all_features_new = col_list_new[:col_list_new.index("Seed_match_interactions_all")]
    print(len(all_features_new))

    col_list_old = list(old.columns)
    all_features_old = col_list_old[:col_list_old.index("Seed_match_interactions_all")]
    print(len(all_features_old))

    # df = set(all_features_new) - set(all_features_old)
    df = set(all_features_old) - set(all_features_new)

    print(df)
    s = list(df)
    df = pd.DataFrame(s)
    df.to_csv(dataset_comper)


# comper()

# comper_intersting()
#
# o = split_features_to_groups_old()
# print(o)
# n = split_features_to_groups_new()
# print(n)
# print(len(o))
# print(len(n))
comper_pos_neg()