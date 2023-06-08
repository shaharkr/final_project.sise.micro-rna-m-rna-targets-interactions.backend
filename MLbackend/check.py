from consts.global_consts import ROOT_PATH,BIOMART_PATH, MERGE_DATA, NEGATIVE_DATA_PATH, GENERATE_DATA_PATH
from utils.utilsfile import read_csv, to_csv
import pandas as pd
from pathlib import Path
from utils.utilsfile import apply_in_chunks, get_wrapper, read_csv, to_csv

def cc():
    file_name = BIOMART_PATH / "human_3utr.csv"
    df_human_3utr = pd.read_csv(file_name)
    df_human_3utr['len'] = df_human_3utr['sequence'].apply(lambda x: len(x))
    df_new = df_human_3utr[df_human_3utr['len']<40]
    print("dd")



def checkMockMirna():
    file_name = ROOT_PATH / "data/positive_interactions/positive_interactions_new/featuers_step/"
    file_name = ROOT_PATH / "data/positive_interactions/positive_interactions_merge/"

    tmp_base = ROOT_PATH / "generate_interactions/mockMirna/"
    features_mock = ROOT_PATH / "data/negative_interactions/mockMirna"
    files = list(features_mock.glob('*features*.csv'))
    run_list = []
    for p in files:
        name = p.name.split('_negative_features.csv')[0] + '_features.csv'
        name = name.split('mockMirna_')[1]
        mock_file = name.split('.csv')[0] + '_negative.csv'
        mock_file = tmp_base / mock_file
        negative_features = features_mock / p
        pos = file_name / name

        duplex_name = p.name.split('_features.csv')[0] + '_duplex.csv'
        duplex_file = features_mock / duplex_name

        normalization_name = p.name.split('_features.csv')[0] + '_normalization.csv'
        normalization_file = features_mock / normalization_name

        df_pos = read_csv(pos)
        # df_pos = df_pos[df_pos['sequence'].apply(lambda x: len(x))<40]
        print(df_pos.shape)

        df_neg_features = read_csv(negative_features)
        df_mock = read_csv(mock_file)
        # d = df_mock['canonic_seed']
        # print(d)
        df_neg_duplex = read_csv(duplex_file)
        df_neg_normalization = read_csv(normalization_file)

        print(df_pos.shape[0] == df_mock.shape[0])
        print(df_pos.shape[0] == df_neg_features.shape[0])
        print("neg_mock_features", df_neg_features.shape[0])
        print("pos", df_pos.shape[0])
        print("mock", df_mock.shape[0])
        print("duplex", df_neg_duplex.shape[0])
        print("normalization", df_neg_normalization.shape[0])
# checkMockMirna()
# s_pos= read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
# print(s_pos.shape)
# s_pos = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/qclash_melanoma_human_ViennaDuplex_features.csv")
# s_neg =read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/mockMirna/mockMirna_qclash_melanoma_human_ViennaDuplex_negative_features.csv")
# # print("s_pos: ", s_pos.shape)
# s_neg = read_csv("/sise/home/efrco/efrco-master/generate_interactions/mockMirna/qclash_melanoma_human_ViennaDuplex_features_negative.csv")
# print("s_neg: ", s_neg.shape)



def checkTarBase():
    tmp_base = ROOT_PATH / "generate_interactions/tarBase/tarBase_human_negative.csv"
    name = "tarBase_human_negative"
    features = ROOT_PATH / "data/negative_interactions/tarBase"
    files = list(features.glob('*features*.csv'))
    for p in files:
        negative_features = features / p
        pos = tmp_base

        duplex_name = name + '_duplex.csv'
        duplex_file = features / duplex_name
        duplex = read_csv(duplex_file)
        normal_name = name + '_normalization.csv'
        normal_file= features / normal_name
        normal = read_csv(normal_file)

        df_pos = read_csv(pos)
        df_neg_features = read_csv(negative_features)
        print(df_pos.shape[0] == df_neg_features.shape[0])
        print("neg_mock_features", df_neg_features.shape[0])
        print("pos", df_pos.shape[0])
        print("duplex:", duplex.shape[0])
        print("normal:", normal.shape[0])



def cleanMock():

    tmp_base = ROOT_PATH / "generate_interactions/mockMirna/"
    features_mock = ROOT_PATH / "data/negative_interactions/mockMirna"
    files = list(features_mock.glob('*features*.csv'))
    run_list = []
    for p in files:
        name = p.name.split('_negative_features.csv')[0] + '_features.csv'
        name = name.split('mockMirna_')[1]
        # mock file
        mock_file = name.split('.csv')[0] + '_negative.csv'
        mock_file = tmp_base / mock_file
        df_mock = read_csv(mock_file)
        df_mock.drop(labels=['canonic_seed', 'non_canonic_seed'], axis=1, inplace=True)
        to_csv(df_mock,mock_file)


        negative_features = features_mock / p
        df_neg_features = read_csv(negative_features)
        df_neg_features.drop(labels=['canonic_seed', 'non_canonic_seed'], axis=1, inplace=True)
        to_csv(df_neg_features, negative_features)


        duplex_name = p.name.split('_features.csv')[0] + '_duplex.csv'
        duplex_file = features_mock / duplex_name
        df_neg_duplex = read_csv(duplex_file)
        df_neg_duplex.drop(labels=['canonic_seed', 'non_canonic_seed'], axis=1, inplace=True)
        to_csv(df_neg_duplex, duplex_file)

        normalization_name = p.name.split('_features.csv')[0] + '_normalization.csv'
        normalization_file = features_mock / normalization_name
        df_neg_normalization = read_csv(normalization_file)
        df_neg_normalization.drop(labels=['canonic_seed', 'non_canonic_seed'], axis=1, inplace=True)
        to_csv(df_neg_normalization, normalization_file)

        # site_name = mockMrna.name.split('_features.csv')[0] + '_site.csv'
        # site_file = features_mock / site_name
        # df_neg_site = read_csv(site_file)
        # df_neg_site.drop(labels=['canonic_seed', 'non_canonic_seed'], axis=1, inplace=True)
        # to_csv(df_neg_site, site_file)



def filter_clash_interaction():

    files_name = ["qclash_melanoma_human_ViennaDuplex_features.csv"]
    df_name = "/sise/home/efrco/efrco-master/generate_interactions/tarBase/tarBase_human_negative.csv"
    df = read_csv(df_name)
    for file_name in files_name:
        file_name_path = MERGE_DATA / file_name
        usecols = ['miRNA ID', 'Gene_ID']
        df_positive_interaction = pd.read_csv(file_name_path, usecols=usecols)

        # transform the format of geneName
        # df_positive_interaction['Gene_ID'] = df_positive_interaction['Gene_ID'].apply(lambda x: x.split("|")[0])
        # df['Gene_ID'] = df['Gene_ID'].apply(lambda x: x.split("|")[0])

        # Intersection to find negative interaction wiche exists in clash poitive interacitons
        intersected_df = pd.merge(df, df_positive_interaction, how='inner', on=['miRNA ID', 'Gene_ID'])

        # remove interactions that exists in both of the dataset
        print("number of rows to remove:" + str(intersected_df.shape[0]))
        print(list(intersected_df['key']))
        print("number of rows before:" + str(df.shape[0]))
        new = df[(~df.key.isin(intersected_df.key))]
        print("number of rows after:" + str(new.shape[0]))
        df = new
        # to_csv(df, df_name)
    return df


def checkTarBase2():
    tmp_base = ROOT_PATH / "generate_interactions/tarBase/tarBase_human_negative.csv"
    name = "tarBase_human_negative"
    features = ROOT_PATH / "data/negative_interactions/tarBase"
    files = list(features.glob('*features*.csv'))
    for p in files:
        negative_features = features / p
        pos = tmp_base

        duplex_name = name + '_duplex.csv'
        duplex_file = features / duplex_name
        duplex = read_csv(duplex_file)
        normal_name = name + '_normalization.csv'
        normal_file = features / normal_name
        normal = read_csv(normal_file)

        df_pos = read_csv(pos)
        df_neg_features = read_csv(negative_features)
        print(df_pos.shape[0] == df_neg_features.shape[0])
        print("neg_mock_features", df_neg_features.shape[0])
        print("pos", df_pos.shape[0])
        print("duplex:", duplex.shape[0])
        print("normal:", normal.shape[0])


#
# file= "/sise/home/efrco/efrco-master/data/negative_interactions/mockMirna/mockMirna_unambiguous_human_ViennaDuplex_negative_normalization.csv"
# file_open= read_csv(file)
# print("Ddddd")


import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def distrubution():


    #################### TarBase ##################################

    file_train_non= Path("/sise/home/efrco/efrco-master/data/train/underSampling/non_overlapping_sites_darnell_human_ViennaDuplex_negative_features_train_underSampling_method.csv")
    df = read_csv(file_train_non)
    # #

    df["len"] = df["site"].apply(lambda x: len(x))
    df["len_seq"] = df["sequence"].apply(lambda x: len(x))
    list_col = ["len"]
    #
    for col in list_col:

        sns.set_theme(style="whitegrid")
        f, ax = plt.subplots(figsize=(9, 9))
        # col_new = df[df[col] < 0.05]
        sns.despine(f)

        sns.histplot(
            df,
            x=col, hue="Label",
            # multiple="stack",
            # palette="light:m_r",
            edgecolor=".10",
            # linewidth=.2,
            # log_scale=False,
        )
        # ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        # ax.set_xticks([0, 0.05, 0.1, 0.15])
        # title = col + "_tarbase"
        # plt.title(title)
        plt.show()
        plt.clf()

    ####################Darnell###############

    # file_train_darnell = Path("/sise/home/efrco/efrco-master/data/train/underSampling/mockMirna_darnell_human_ViennaDuplex_negative_features_train_underSampling_method.csv")
    # df = read_csv(file_train_darnell)


    # for col in list_col:
    #     # col_new = df[df[col] < 0.05]
    #     sns.set_theme(style="whitegrid")
    #     f, ax = plt.subplots(figsize=(9, 9))
    #     sns.despine(f)
    #
    #     sns.histplot(
    #         df,
    #         x=col, hue="Label",
    #         # multiple="stack",
    #         # palette="light:m_r",
    #         edgecolor=".10",
    #         # linewidth=.2,
    #         # log_scale=False,
    #     )
    #     # ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    #     # ax.set_xticks([0, 0.05, 0.1, 0.15])
    #     title = col + "_darnell"
    #     plt.title(title)
    #     plt.show()
    #     plt.savefig(f'/sise/home/efrco/efrco-master/data/fig/{col}_darnell.png', format="png")
    #
    #     plt.clf()

# distrubution()


def train_test_balance():

    file_train_tarBase = Path("/sise/home/efrco/efrco-master/data/train/underSampling/tarBase_human_negative_features_train_underSampling_method.csv")
    df_tarBase = read_csv(file_train_tarBase)

    file_train_darnell = Path("/sise/home/efrco/efrco-master/data/train/underSampling/mockMirna_darnell_human_ViennaDuplex_negative_features_train_underSampling_method.csv")
    df_darnell = read_csv(file_train_darnell)

    df_tarBase_pos = df_tarBase[df_tarBase['Label'] == 1]
    print(df_tarBase_pos.shape)
    df_darnell_pos = df_darnell[df_darnell['Label'] == 1]
    new=df_darnell_pos == df_tarBase_pos
    print(df_darnell_pos == df_tarBase_pos)


def mirna_dist():
    path_tarBase = "/sise/home/efrco/efrco-master/data/negative_interactions/tarBase/tarBase_human_negative_features.csv"
    df_tarBase = read_csv(path_tarBase)
    count_mirna = 0

    path_darnell = "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv"
    df_darnell = read_csv(path_darnell)


    mirna_tarBase = df_tarBase.groupby(['miRNA ID']).size().reset_index(name='counts_tarBase')
    count_mirna = count_mirna - mirna_tarBase.shape[0] #240

    # print(mirna_tarBase)

    mirna_darnell = df_darnell.groupby(['miRNA ID']).size().reset_index(name='counts_darnell')
    count_mirna = count_mirna - mirna_darnell.shape[0] #222


    # print(mirna_darnell)
    dist = pd.merge(mirna_tarBase, mirna_darnell, how='inner', on=['miRNA ID'])
    # dist.rename(columns={"counts_x": "counts_tarBase", "counts_y": "counts_darnell"}, inplace= True)
    count_mirna = count_mirna + dist.shape[0]  #95

    mirna_list_tarBase = list(mirna_tarBase['miRNA ID'])
    mirna_list_darnell = list(mirna_darnell['miRNA ID'])

    for mir in mirna_list_tarBase:
        if mir not in mirna_list_darnell:
            count_mirna += 1
            row = mirna_tarBase[mirna_tarBase['miRNA ID'] == mir]
            cont_mir = row.iloc[0]['counts_tarBase']

            dist = dist.append({'miRNA ID': mir, "counts_tarBase": cont_mir, "counts_darnell": 0}, ignore_index=True)

    for mir in mirna_list_darnell:
        if mir not in mirna_list_tarBase:
            count_mirna += 1
            row = mirna_darnell.loc[mirna_darnell['miRNA ID'] == mir]
            cont_mir = row.iloc[0]['counts_darnell']
            dist = dist.append({'miRNA ID': mir, "counts_tarBase": 0,"counts_darnell": cont_mir}, ignore_index=True)

    print(count_mirna)
    print(sum(dist['counts_darnell']) == df_darnell.shape[0])
    print(sum(dist['counts_tarBase']) == df_tarBase.shape[0])


    full_path = ROOT_PATH / "distribution_mirna_new.csv"
    to_csv(dist, full_path)


def mirna_dist_seed():
    path_tarBase = "/sise/home/efrco/efrco-master/data/negative_interactions/tarBase/tarBase_human_negative_features.csv"
    df_tarBase = read_csv(path_tarBase)

    path_darnell = "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/darnell_human_ViennaDuplex_features.csv"
    df_darnell = read_csv(path_darnell)


    mirna_tarBase = df_tarBase.groupby(['seed_family']).size().reset_index(name='counts_tarBase')

    # print(mirna_tarBase)

    mirna_darnell = df_darnell.groupby(['seed_family']).size().reset_index(name='counts_darnell')


    # print(mirna_darnell)
    dist = pd.merge(mirna_tarBase, mirna_darnell, how='inner', on=['seed_family'])
    # dist.rename(columns={"counts_x": "counts_tarBase", "counts_y": "counts_darnell"}, inplace= True)

    seed_list_tarBase = list(mirna_tarBase['seed_family'])
    seed_list_darnell = list(mirna_darnell['seed_family'])

    for seed in seed_list_tarBase:
        if seed not in seed_list_darnell:
            row = mirna_tarBase[mirna_tarBase['seed_family'] == seed]
            cont_seed = row.iloc[0]['counts_tarBase']

            dist = dist.append({'seed_family': seed, "counts_tarBase": cont_seed, "counts_darnell": 0}, ignore_index=True)

    for seed in seed_list_darnell:
        if seed not in seed_list_tarBase:
            row = mirna_darnell.loc[mirna_darnell['seed_family'] == seed]
            cont_seed = row.iloc[0]['counts_darnell']
            dist = dist.append({'seed_family': seed, "counts_tarBase": 0,"counts_darnell": cont_seed}, ignore_index=True)

    print(sum(dist['counts_darnell']) == df_darnell.shape[0])
    print(sum(dist['counts_tarBase']) == df_tarBase.shape[0])


    full_path = ROOT_PATH / "distribution_mirna_seed_new.csv"
    to_csv(dist, full_path)




# train_test_balance()
# mirna_dist()
# mirna_dist_seed()


def size_mock():

    file1_before = "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/data_without_featuers/darnell_human_ViennaDuplex.csv"
    file1_after = "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv"
    f1 = read_csv(file1_before)
    f2 = read_csv(file1_after)
    f3= read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/darnell_human_ViennaDuplex_features.csv")
    print(f1.shape)
    print(f2.shape)
    print(f3.shape)
#
# size_mock()

#
# file = read_csv(Path("/sise/home/efrco/efrco-master/data/negative_interactions/non_overlapping_sites/non_overlapping_sites_darnell_human_ViennaDuplex_negative_features.csv"))
# file = read_csv(Path("/sise/home/efrco/efrco-master/data/negative_interactions/non_overlapping_sites/non_overlapping_sites_darnell_human_ViennaDuplex_negative_features.csv"))
# file["len"] = file["site"].apply(lambda x: len(x))
# print(file.shape)

# file["len_seq"] = file["full_mrna"].apply(lambda x: len(x))

# print("f")
# sns.set_theme(style="whitegrid")
# f, ax = plt.subplots(figsize=(9, 9))
# # col_new = df[df[col] < 0.05]
# sns.despine(f)
#
# sns.histplot(
#     file,
#     x="len",
#     # multiple="stack",
#     # palette="light:m_r",
#     edgecolor=".10",
#     # linewidth=.2,
#     # log_scale=False,
# )
# ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
# ax.set_xticks([0, 0.05, 0.1, 0.15])
# plt.show()
# plt.clf()
#
#
#
# print(file.shape)



def checkMockMrna():
    file_name = ROOT_PATH / "data/positive_interactions/positive_interactions_new/featuers_step/"
    tmp_base = ROOT_PATH / "generate_interactions/mockMrna/"
    features_mock = ROOT_PATH / "data/negative_interactions/mockMirna"
    files = list(features_mock.glob('*features*.csv'))
    run_list = []
    for p in files:
        name = p.name.split('_negative_features.csv')[0] + '_features.csv'
        name = name.split('mockMirna_')[1]
        mock_file = name.split('.csv')[0] + '_negative.csv'
        mock_file = tmp_base / mock_file
        negative_features = features_mock / p
        pos = file_name / name

        duplex_name = p.name.split('_features.csv')[0] + '_duplex.csv'
        duplex_file = features_mock / duplex_name

        normalization_name = p.name.split('_features.csv')[0] + '_normalization.csv'
        normalization_file = features_mock / normalization_name

        df_pos = read_csv(pos)
        # df_pos = df_pos[df_pos['sequence'].apply(lambda x: len(x))<40]
        print(df_pos.shape)

        df_neg_features = read_csv(negative_features)
        df_mock = read_csv(mock_file)
        # d = df_mock['canonic_seed']
        # print(d)
        df_neg_duplex = read_csv(duplex_file)
        df_neg_normalization = read_csv(normalization_file)

        print(df_pos.shape[0] == df_neg_features.shape[0])
        print("neg_mock_features", df_neg_features.shape[0])
        print("pos", df_pos.shape[0])
        print("mock", df_mock.shape[0])
        print("duplex", df_neg_duplex.shape[0])
        print("normalization", df_neg_normalization.shape[0])
# checkMockMrna()
# s_pos= read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
# print(s_pos.shape)
# s_pos = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/qclash_melanoma_human_ViennaDuplex_features.csv")
# s_neg =read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/mockMirna/mockMirna_qclash_melanoma_human_ViennaDuplex_negative_features.csv")
# # print("s_pos: ", s_pos.shape)
# s_neg = read_csv("/sise/home/efrco/efrco-master/generate_interactions/mockMirna/qclash_melanoma_human_ViennaDuplex_features_negative.csv")
# print("s_neg: ", s_neg.shape)

# def insert_mock_site(full_mrna, start, end, site):
#     full_mrna = full_mrna[:start] + site + full_mrna[end + 1:]
#
#     return full_mrna
#
#
# print(insert_mock_site("efrat",1,2,"mmm"))

def null_in_dataset_h3():
    my = read_csv(
        "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
    C1 = my[my.isna().any(axis=1)]
    my = read_csv(
        "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
    C1 = my[my.isna().any(axis=1)]
    my_after_drop_na = my.dropna()

    # GILAD
    original = read_csv(
        "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/darnell_human_ViennaDuplex_features.csv")
    C2 = original[original.isna().any(axis=1)]
    # GILAD
    original = read_csv(
        "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/darnell_human_ViennaDuplex_features.csv")
    original.drop(columns=['miRNAMatchPosition_21','miRNAMatchPosition_22'], inplace=True)
    C2 = original[original.isna().any(axis=1)]
    C2['size'] = C2['miRNA sequence'].apply(lambda x: len(x))
    print(C1.shape)



def null_in_dataset_tarBase():
    my = read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/tarBase/tarBase_human_negative_features.csv")
    C1 = my[my.isna().any(axis=1)]
    my_after_drop_na = my.dropna()
    cols = my.columns[my.isna().any()].tolist()
    print(cols)
    path = "/sise/home/efrco/efrco-master/tarbasenull.csv"
    to_csv(C1, path)



# null_in_dataset_tarBase()




# file_name = ROOT_PATH / "data/positive_interactions/positive_interactions_merge"
# files = list(file_name.glob('**/*.csv'))
# for p in files:
#     df = read_csv(p)
#     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#     print(p)
#     print(df.shape)
#     C1 = df[df.isna().any(axis=1)]
#     print(C1.shape)

#
# import Classifier.FeatureReader as FeatureReader
# from Classifier.FeatureReader import get_reader
#
# FeatureReader.reader_selection_parameter = "without_hot_encoding"
#
# feature_reader = get_reader()
#
# X, y = feature_reader.file_reader("/sise/home/efrco/efrco-master/data/test/one_class_svm/non_overlapping_sites_darnell_human_ViennaDuplex_negativeone_class.csv")
# print("f")

# Example of a confusion matrix in Python
# from sklearn.metrics import confusion_matrix
#
# expected = [1,0,1,1]
# predicted = [0,1,1,1]
# results = confusion_matrix(expected, predicted)
# print(results.ravel())
# print("(tn, fp, fn, tp)")
# # print(results)
# TP = results[0][0]
# FP = results[0][1]
# FN = results[1][0]
# TN = results[1][1]
# print("TP:", TP)
# print("FP:", FP)
# print("FN:", FN)
# print("TN:", TN)
# import seaborn as sns
#
# ax = sns.heatmap(results, annot=True, cmap='Blues')
#
# ax.set_title('Seaborn Confusion Matrix with labels\n\n');
# ax.set_xlabel('\nPredicted Values')
# ax.set_ylabel('Actual Values ');
#
# ## Ticket labels - List must be in alphabetical order
# ax.xaxis.set_ticklabels(['False','True'])
# ax.yaxis.set_ticklabels(['False','True'])
#
# ## Display the visualization of the Confusion Matrix.
# plt.show()


# path ="/sise/home/efrco/efrco-master/data/negative_interactions/mockMrna/mockMrna_nucleotides_method1_darnell_human_ViennaDuplex_negative_features.csv"
path = "/sise/home/efrco/efrco-master/data/negative_interactions/mockMrna/mockMrna_denucleotides_method3_darnell_human_ViennaDuplex_negative_features.csv"
path = "/sise/home/efrco/efrco-master/data/negative_interactions/mockMirna/mockMirna_darnell_human_ViennaDuplex_negative_features.csv"
path = "/sise/home/efrco/efrco-master/data/negative_interactions/non_overlapping_sites/non_overlapping_sites_darnell_human_ViennaDuplex_negative_features.csv"
path = "/sise/home/efrco/efrco-master/data/negative_interactions/tarBase/tarBase_human_negative_features.csv"
# neg = read_csv(path)
# neg = neg[neg.isna().any(axis=1)]
# cols = neg.columns[neg.isna().any()].tolist()
# print(cols)
# print(neg.shape)
# path_output = "/sise/home/efrco/efrco-master/denucleotides_method2_null_neg.csv"
# to_csv(neg, path_output)


def creat_dir():
    import os

    # Directory
    for i in range(12,13):
        directory = "split_file"
        name = "clip_" + str(i)
        # Parent Directory path
        parent_dir = Path("/sise/vaksler-group/IsanaRNA/miRNA_target_rules/Isana/") / name
        # parent_dir = "/sise/home/efrco/efrco-master/Results/models/models_binary_comper/"


        # Path
        path = os.path.join(parent_dir, directory)

        # Create the directory
        # 'GeeksForGeeks' in
        # '/home / User / Documents'
        os.mkdir(path)
        print("Directory '% s' created" % directory)

# creat_dir()

# def model_measuremnt_try():
#     from sklearn.datasets import make_classification
#     from sklearn.linear_model import LogisticRegression
#     from sklearn.model_selection import train_test_split
#     from sklearn.metrics import precision_recall_curve
#     from sklearn.metrics import f1_score
#     from sklearn.metrics import auc
#     from matplotlib import pyplot
#     # generate 2 class dataset
#     X, y = make_classification(n_samples=1000, n_classes=2, random_state=1)
#     # split into train/test sets
#     trainX, testX, trainy, testy = train_test_split(X, y, test_size=0.5, random_state=2)
#     # fit a model
#     model = LogisticRegression(solver='lbfgs')
#     model.fit(trainX, trainy)
#     # predict probabilities
#     lr_probs = model.predict_proba(testX)
#     # keep probabilities for the positive outcome only
#     lr_probs = lr_probs[:, 1]
#     # predict class values
#     yhat = model.predict(testX)
#     lr_precision, lr_recall, _ = precision_recall_curve(testy, lr_probs)
#     lr_f1, lr_auc = f1_score(testy, yhat), auc(lr_recall, lr_precision)
#     # summarize scores
#     print('Logistic: f1=%.3f auc=%.3f' % (lr_f1, lr_auc))
#     # plot the precision-recall curves
#     no_skill = len(testy[testy == 1]) / len(testy)
#     pyplot.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
#     pyplot.plot(lr_recall, lr_precision, marker='.', label='Logistic')
#     # axis labels
#     pyplot.xlabel('Recall')
#     pyplot.ylabel('Precision')
#     # show the legend
#     pyplot.legend()
#     # show the plot
#     pyplot.show()
#
# model_measuremnt_try()

# mirna = read_csv("/sise/vaksler-group/IsanaRNA/miRNA_target_rules/Isana/clip_3/mirna.csv")



def generate_up(number):
    list_numbrers = []
    for i in range(number,number+5):
        list_numbrers.append(i)
    return list_numbrers


def generate_down(number):
    list_numbrers = []
    for i in range(number-4, number+1):
        list_numbrers.append(i)
    return list_numbrers



def get_max_str(lst):
    # in case that we don't found overlap between two sequences
    if lst == []:
        return ""
    return max(lst, key=len)

def substringFinder(string1, string2):
    answer = ""
    anslist = []
    len1, len2 = len(string1), len(string2)
    for i in range(len1):
        match = ""
        for j in range(len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                #if (len(match) > len(answer)):
                answer = match
                if answer != '' and len(answer) > 1:
                    anslist.append(answer)
                match = ""

        if match != '':
            anslist.append(match)
    c = get_max_str(anslist)
    return get_max_str(anslist)





# df = mrna.loc[df_g["sequence length"].idxmax()]

# print("number:", df_g.ngroups)
# print(df.shape)
# # to_csv(df, path)
# #
#
# path = "/sise/home/efrco/efrco-master/generate_interactions/clip_interaction/clip_3.csv"
# df = read_csv(path)
# print(df.shape)

from consts.global_consts import ROOT_PATH, DATA_PATH,CLIP_PATH_DATA

def mirna_dist_clip_interactions():
    clip_data_path = CLIP_PATH_DATA
    df = []
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*mrna.csv*"):
            mirna_df = read_csv(clip_dir / "mirna.csv")
            df.append(mirna_df)
    result = pd.concat(df)
    result = result.groupby(['miRNA ID']).size().reset_index(name='count')
    # result = result.groupby(['miRNA ID','count']).size().reset_index(name='count of mirna')
    # print(result)
    # # result = result[result['count'] > 12]
    # ax = sns.barplot(data=result, x="count", y="count of mirna")
    # ax.bar_label(ax.containers[0])
    # plt.show()
    return result


def mrna_dist_clip_interactions():
    clip_data_path = CLIP_PATH_DATA
    df = []
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*mrna_clean.csv*"):
            mirna_df = read_csv(clip_dir / "mrna_clean.csv")
            df.append(mirna_df)
    result = pd.concat(df)
    result = result.groupby(['ID']).size().reset_index(name='count')
    # result = result.groupby(['count']).size().reset_index(name='count of mrna')
    print(result)
    # result = result[result['count'] > 20]
    ax = sns.barplot(data=result, x="ID", y="count")
    ax.bar_label(ax.containers[0])
    path= "/sise/home/efrco/mrna_dist.csv"
    plt.show()
    to_csv(result, path)

# mrna_dist_clip_interactions()
# mirna_dist_clip_interactions()

def clip_files_dist_mirna():
    clip3_df = read_csv("/sise/home/efrco/efrco-master/generate_interactions/clip_interaction/clip_6.csv")
    print("f")
    result = mirna_dist_clip_interactions()
    result = result[result['count'] == 2]
    mirna_list = result['miRNA ID'].tolist()
    print(len(mirna_list))
    clip3_df['exists'] = clip3_df['miRNA ID'].apply(lambda id: 1 if id in mirna_list else 0)
    mirna_list_clip = clip3_df['miRNA ID'].tolist()
    print(clip3_df['exists'].value_counts())


    count = 0
    for mir in mirna_list:
        if mir not in mirna_list_clip:
            print("BUG")
        else:
            count = count + 1
    print(count)

def dist_enrgy_positive_negative():

    # file_pos = read_csv(Path("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv"))
    # file_neg = read_csv(Path("/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_features.csv"))
    r = read_csv("/sise/home/efrco/efrco-master/data/train/underSampling/0/clip_interaction_clip_3_negative_features_train_underSampling_method_0.csv")
    file_pos = r[r['Label'] == 1]
    file_neg = r[r['Label'] == 0]
    sns.set_theme(style="whitegrid")
    f, ax = plt.subplots(figsize=(9, 9))
    # col_new = df[df[col] < 0.05]
    sns.despine(f)

    sns.histplot(
        file_pos,
        x="Energy_MEF_Duplex",
    )
    # ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.set_xticks([0, 0.05, 0.1, 0.15])
    plt.show()
    plt.clf()
    # print(file.shape)

# dist_enrgy_positive_negative()

def dist_enrgy_positive_negative_same_graph():

    r_clip = read_csv("/sise/home/efrco/efrco-master/data/train/underSampling/0/clip_interaction_clip_3_negative_features_train_underSampling_method_0.csv")
    file_neg_clip = r_clip[r_clip['Label'] == 0]
    file_neg_clip['method'] = 'neg_clip'
    file_neg_clip.reset_index(drop=True, inplace=True)
    file_neg_clip["len"] = file_neg_clip["site"].apply(lambda x: len(x))


    r_mock = read_csv("/sise/home/efrco/efrco-master/data/train/underSampling/0/mockMirna_darnell_human_ViennaDuplex_negative_features_train_underSampling_method_0.csv")
    file_neg_mock = r_mock[r_mock['Label'] == 0]
    file_neg_mock['method'] = 'mock'
    file_neg_mock.reset_index(drop=True, inplace=True)

    file_pos = r_clip[r_clip['Label'] == 1]
    file_pos["len"] = file_pos["site"].apply(lambda x: len(x))
    file_neg_mock["len"] = file_neg_mock["site"].apply(lambda x: len(x))
    file_pos['method'] = 'pos'

    merge_data_frame= pd.concat([file_neg_mock, file_neg_clip, file_pos])
    print("pos", merge_data_frame[merge_data_frame['Label']==1].shape)
    print("neg", merge_data_frame[merge_data_frame['Label']==0].shape)

    merge_data_frame.reset_index(inplace=True)
    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(7, 5))
    sns.despine(f)

    # sns.histplot(
    #     merge_data_frame,
    #     x="Energy_MEF_Duplex", hue="method",
    #     multiple="stack",
    #     palette="light:m_r",
    #     edgecolor=".3",
    #     stat='density'
    #     # linewidth=.5,
    #     # log_scale=True,
    # )
    sns.ecdfplot(data=merge_data_frame, x="Energy_MEF_Duplex", hue="method")

    plt.show()
    plt.clf()

# dist_enrgy_positive_negative_same_graph()


def dist_enrgy_positive_negative_same_graph_len_site_original():

    pos = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/data_without_featuers/darnell_human_ViennaDuplex.csv")
    pos['method'] = 'pos'
    pos["len"] = pos["site"].apply(lambda x: len(x))
    pos.reset_index(drop=True, inplace=True)

    neg = read_csv("/sise/home/efrco/efrco-master/generate_interactions/clip_interaction/clip_3.csv")
    neg["len"] = neg["site"].apply(lambda x: len(x))
    neg['method'] = 'neg'
    neg.reset_index(drop=True, inplace=True)


    merge_data_frame= pd.concat([pos, neg])

    merge_data_frame.reset_index(inplace=True)
    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(7, 5))
    sns.despine(f)
    sns.ecdfplot(data=merge_data_frame, x="len", hue="method")

    plt.show()
    plt.clf()

# dist_enrgy_positive_negative_same_graph_len_site_original()

def dist_enrgy_positive_negative_same_graph_len_site_of_the_new_site_after_duplex():

    pos =read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex.csv")
    pos['method'] = 'pos'
    pos["len"] = pos["site"].apply(lambda x: len(x))
    pos.reset_index(drop=True, inplace=True)

    neg = read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_duplex.csv")
    neg["len"] = neg["site"].apply(lambda x: len(x))
    neg['method'] = 'neg'
    neg.reset_index(drop=True, inplace=True)


    merge_data_frame= pd.concat([pos, neg])

    merge_data_frame.reset_index(inplace=True)
    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(7, 5))
    sns.despine(f)
    sns.ecdfplot(data=merge_data_frame, x="len", hue="method")

    plt.show()
    plt.clf()

# dist_enrgy_positive_negative_same_graph_len_site_of_the_new_site_after_duplex()


def dist_enrgy_positive_negative_same_graph_number_mis_tail_site():

    pos_duplex = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/duplex_step/darnell_human_ViennaDuplex_duplex.csv")
    pos_duplex['method'] = 'pos'

    neg_duplex_clip= read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_duplex.csv")
    neg_duplex_clip['method'] = 'neg_duplex_clip_without_extend'

    neg_duplex_to_original_site= read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_d_original.csv")
    neg_duplex_to_original_site['method'] = 'neg_duplex_clip_with_extend'


    merge_data_frame= pd.concat([pos_duplex, neg_duplex_clip,neg_duplex_to_original_site])

    merge_data_frame.reset_index(inplace=True)

    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(7, 5))
    sns.despine(f)
    # sns.ecdfplot(data=merge_data_frame, x="not_match_site", hue="method")
    plt.show()
    plt.clf()

# dist_enrgy_positive_negative_same_graph_number_mis_tail_site()


def expermient():

    df_clip = read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_features.csv")
    df_clip = df_clip[df_clip['not_match_site'] < 5]
    df_clip = df_clip[df_clip['len_site'] > 74]
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(df_clip.shape)

    path_neg = "/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_features.csv"
    to_csv(df_clip, path_neg)


    df_pos = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex.csv")
    df_pos = df_pos[df_pos['not_match_site'] < 5]
    path_pos = "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv"
    to_csv(df_pos, path_pos)

# expermient()

def add_not_match_site_to_non_extend_site_clip_3():
    fin = "/sise/home/efrco/efrco-master/generate_interactions/clip_interaction/clip_3.csv"
    fout = "/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_duplex.csv"
    from pipeline_steps.duplex_step import duplex as duplex_positive
    # negative interactions
    duplex_positive('ViennaDuplex', fin, fout)
# add_not_match_site_to_non_extend_site_clip_3()

from consts.global_consts import MERGE_DATA, ROOT_PATH, BIOMART_PATH, POSITIVE_PATH_FEATUERS

def isana_task(name_gene):

    files_name = ["human_mapping_ViennaDuplex_features.csv", "darnell_human_ViennaDuplex_features.csv","unambiguous_human_ViennaDuplex_features.csv", "qclash_melanoma_human_ViennaDuplex_features.csv"]

    for file_name in files_name:
        file_name_path = POSITIVE_PATH_FEATUERS / file_name
        df_positive_interaction = pd.read_csv(file_name_path)
        df_positive_interaction['len'] = df_positive_interaction['site'].apply(lambda x: len(x))

        # transform the format of geneName
        df_positive_interaction['Gene_ID'] = df_positive_interaction['Gene_ID'].apply(lambda x: x.split("|")[0])
        df_positive_filter = df_positive_interaction[df_positive_interaction['Gene_ID'] == name_gene]
        df_positive_filter.drop(columns= ['level_0'],inplace=True)
        path = Path("/sise/home/efrco/")
        name_file_new = name_gene + "_" + file_name
        path_write = path / name_file_new
        # to_csv(df_positive_filter, path_write )


# isana_task(name_gene="ENSG00000079999")




def read(name_gene):
    # Consts
    file_name = ROOT_PATH / "generate_interactions/tarBase/interaction.xlsx"
    validation_string = "geneID"
    # usecols = ['geneID', 'geneName',
    #            'mirna', 'species', 'tissue', 'method']

    df = pd.read_excel(file_name, engine='openpyxl')
    assert df.columns[0] == validation_string, f"reader validation error: {df.columns[0]}"
    print("Number of interactions in tarBase is:", df.shape[0])

    # # filter only liver interaction that will be match to darnell
    df_filter = df[df['positive_negative'] == 'POSITIVE']
    df_filter = df_filter[df_filter['geneID'] ==name_gene]
    df_filter = df_filter[df_filter['species'] == 'Homo sapiens']

    path = Path("/sise/home/efrco/")
    name_file_new = name_gene + "_" + "tarBase.csv"

    path_write = path / name_file_new
    to_csv(df_filter, path_write)
    # print("Number of after filter is:", df.shape[0])

    return df


def tail_check():
    clash_darnell = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
    clip = read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/clip_interaction/clip_interaction_clip_3_negative_features.csv")
    print("f")

# tail_check()

import numpy as np


def check_df_non(Seed_match_noncanonical, non_canonic_seed):
    if non_canonic_seed == Seed_match_noncanonical:
        return True
    else:
        return False


def check_df_canon(non_canonic_seed, Seed_match_canonical):
    if (non_canonic_seed == Seed_match_canonical):
        return True
    else:
        return False
def canon_noncanon():

    usecols_gilad = ['mRNA_name', 'target sequence','microRNA_name', 'non_canonic_seed', 'canonic_seed']
    d = read_csv("/sise/vaksler-group/IsanaRNA/miRNA_target_rules/benorgi/TPVOD/Data/Features/CSV/human_dataset3_duplex_positive_feature.csv")
    gilad_data = pd.read_csv("/sise/vaksler-group/IsanaRNA/miRNA_target_rules/benorgi/TPVOD/Data/Features/CSV/human_dataset3_duplex_positive_feature.csv",  usecols=usecols_gilad)
    gilad_data.rename(columns={"microRNA_name": "miRNA ID", "target sequence":"site"}, inplace=True)
    gilad_data['site'] = gilad_data['site'].replace(to_replace='T', value='U', regex=True)
    gilad_data['mRNA_name'] = gilad_data['mRNA_name'].apply(lambda x: x.split("|")[0:2])
    gilad_data['ID'] = gilad_data['mRNA_name'].apply(lambda x: x[0])
    gilad_data['mRNA_name'] = gilad_data['mRNA_name'].apply(lambda x: x[0] + "|" + x[1])

    gilad_data.rename(columns={"mRNA_name": "Gene_ID"}, inplace=True)


    usecols_efrat = ['Gene_ID', 'site','miRNA ID', 'Seed_match_noncanonical', 'Seed_match_canonical']
    clash_darnell = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
    efrat_data = pd.read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/darnell_human_ViennaDuplex_features.csv",  usecols=usecols_efrat)
    efrat_data['ID'] = efrat_data['Gene_ID'].apply(lambda x: x.split("|"))
    efrat_data['ID'] = efrat_data['ID'].apply(lambda x: x[0])


    intersected_df = pd.merge(gilad_data, efrat_data, how='inner', on=['miRNA ID','ID'] )

    intersected_df['non_canonic_result'] = intersected_df.apply(func=get_wrapper(check_df_non,
                                           "Seed_match_noncanonical", "non_canonic_seed"), axis=1)
    intersected_df['canonic_result'] = intersected_df.apply(func=get_wrapper(check_df_canon,
                                           "Seed_match_canonical", "canonic_seed",
                                         ), axis=1)

    problem   =intersected_df[intersected_df['canonic_result'] == False]
    print(problem.shape) #1214 probelm

# canon_noncanon()

