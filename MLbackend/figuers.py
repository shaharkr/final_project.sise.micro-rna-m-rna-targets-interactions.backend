from consts.global_consts import ROOT_PATH,BIOMART_PATH, MERGE_DATA, NEGATIVE_DATA_PATH, GENERATE_DATA_PATH
from utils.utilsfile import read_csv, to_csv
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Classifiers_Cross_Validiation_Runing import CrossValidation

path_figures = Path("/sise/home/efrco/efrco-master/figuers/")


def canon_distribution_figures_dataset(data_set_pos):
    dir = NEGATIVE_DATA_PATH
    all_frames = []
    for method_dir in dir.iterdir():
        print(method_dir)
        list_method = ['mockMirna', 'non_overlapping_sites']

        for dataset_file_neg in method_dir.glob("*features*"):
            if any(method in method_dir.stem for method in list_method):
                correct_dataset = data_set_pos.split("_features")[0] + "_negative_features"
                list_split = dataset_file_neg.stem.split("/")
                list_split = [str(s) for s in list_split]
                if correct_dataset not in list_split[0]:
                    continue
            neg = read_csv(dataset_file_neg)
            neg['dataset_name'] = method_dir.name
            all_frames.append(neg)
            # break
        # break

    # add darnell
    # path_pos = read_csv("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
    path_pos = read_csv("/sise/vaksler-group/IsanaRNA/miRNA_target_rules/benorgi/TPVOD/Data/Features/CSV/human_dataset3_duplex_positive_feature.csv")
    path_pos['dataset_name'] = 'positive_interactions'
    all_frames.append(path_pos)

    result = pd.concat(all_frames)
    result = result[['Seed_match_canonical', 'Seed_match_noncanonical', 'dataset_name']]
    result['type_seed'] = result.apply(lambda row: 'canonical' if row['Seed_match_canonical'] else 'noncanonical', axis = 1)
    result.set_index('dataset_name')
    cross_tab_prop = pd.crosstab(index=result['dataset_name'],
                                 columns=result['type_seed'],
                                 normalize="index")
    cross_tab = pd.crosstab(index=result['dataset_name'],
                            columns=result['type_seed'])

    cross_tab_prop.plot(kind='bar',
                        stacked=True,
                        colormap='tab10',
                        figsize=(10, 6), color=['red', 'skyblue'])

    for n, x in enumerate([*cross_tab.index.values]):
        for (proportion, y_loc) in zip(cross_tab_prop.loc[x],
                                       cross_tab_prop.loc[x].cumsum()):
            plt.text(x=n - 0.17,
                     y=(y_loc - proportion) + (proportion / 2),
                     s=f'{np.round(proportion * 100, 1)}%',
                     color="black",
                     fontsize=12,
                     fontweight="bold")

    for n, x in enumerate([*cross_tab.index.values]):
        table = cross_tab.loc[x]
        canon = table['canonical']
        noncanon = table['noncanonical']
        total = canon + noncanon

        plt.text(x=n - 0.17,
                 y=y_loc,
                 s=f'{total}',
                 color="black",
                 fontsize=12,
                 fontweight="bold")


    plt.legend(loc='lower left', ncol=2)
    plt.xlabel("Dataset")
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Percentage")
    fname = path_figures / "canon_distribution_figures_dataset.pdf"
    plt.savefig(fname, format="pdf", bbox_inches='tight')
    plt.show()
    plt.clf()


# canon_distribution_figures_dataset(data_set_pos="darnell_human_ViennaDuplex_features")


def creat_figures_Intra_anaylsis():
    cv = CrossValidation("darnell_human_ViennaDuplex_features", "measurments_cross_validation", number_iterations=15)
    list_metrix = ['ACC', 'TPR', 'TNR', 'PPV']

    #####################TarBase############################
    # for measuer in list_metrix:
    #     cv.summary_matrix_tarBase(measuer, target_calculate= "std")
    #     cv.summary_matrix_tarBase(measuer, target_calculate= "mean")
    #
    # #####################MockMrna############################
    #
    # for measuer in list_metrix:
    #     cv.summary_matrix_mock_mrna(measuer, target_calculate= "std")
    #     cv.summary_matrix_mock_mrna(measuer, target_calculate= "mean")

    #####################non_overlapping_sites############################

    # for measuer in list_metrix:
    #     cv.summary_matrix_non_overlapping_site(measuer, target_calculate= "std")
    #     cv.summary_matrix_non_overlapping_site(measuer, target_calculate= "mean")

    #####################clip_data_tail############################
    #
    # for measuer in list_metrix:
    #     cv.summary_matrix_clip(measuer, target_calculate="std")
    #     cv.summary_matrix_clip(measuer, target_calculate="mean")

    cv.summary_matrix_Intra("ACC")



# creat_figures_Intra_anaylsis()

# d = read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/non_overlapping_sites/non_overlapping_sites_random_darnell_human_ViennaDuplex_negative_features.csv")
# print("f")

df = read_csv("/sise/home/efrco/efrco-master/data/negative_interactions/tarBase/tarBase_Liver_human_negative_features.csv")
print(df[col].unique())
