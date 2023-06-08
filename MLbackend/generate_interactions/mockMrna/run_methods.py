from pipline_steps_negative.duplex_step_negative import duplex as duplex_negative
from pipline_steps_negative.rna_site_insertion_negative import get_site_from_extended_site
from pipline_steps_negative.normalization_final_step_negative import finalize
from pipeline_steps.feature_extraction import feature_extraction
from consts.global_consts import NEGATIVE_DATA_PATH, GENERATE_DATA_PATH


# step 1- formalization all the dataset for the same format


def full_pipline(name_of_method: str, name_shuffle: str , name_of_file: str, duplex=duplex_negative):
    # step 2- extract duplex of the interaction by VieannaDuplex
    # only for mock mirna
    name_of_file_primary = name_of_method + "_" + name_shuffle + "_" + name_of_file.split('_features')[0] + "_negative"

    # name_of_file_primary = name_of_method + "_" + name_of_file
    fout_primary = NEGATIVE_DATA_PATH / name_of_method

    name_of_file = name_of_file.split('_negative')[0] + "_" + name_shuffle + "_negative.csv"
    fin = GENERATE_DATA_PATH / name_of_method / name_of_file

    name_of_file = name_of_file_primary + "_duplex.csv"
    fout= fout_primary / name_of_file

    # negative interactions
    print("###############Duplex NEGATIVE#############")
    duplex('ViennaDuplex', fin, fout)

    # step 3- extract the site and his coordination's
    fin = fout
    print("###############Site#############")

    name_of_file = name_of_file_primary + "_site.csv"
    fout = fout_primary /  name_of_file
    get_site_from_extended_site(fin, fout)

    print("###############Normaliztion#############")

    # step 4- normalization of the dataframe
    fin = fout
    name_of_file = name_of_file_primary + "_normalization.csv"
    fout = fout_primary / name_of_file
    finalize(fin, fout)

    # step 5- extract features
    print("###############extract features#############")

    fin = fout
    name_of_file = name_of_file_primary + "_features.csv"
    fout = fout_primary / name_of_file
    name_of_file_ch = name_of_file_primary + "_check_before_filter.csv"
    fout_check = fout_primary / name_of_file_ch
    feature_extraction(fin, fout)

def run():
    ####### mockMrna 1 approach ######
    #5266683

    # from generate_interactions.mockMrna import NegativeMockMrna
    # NegativeMockMrna.main(name_shuffle='nucleotides')
    # full_pipline(name_of_method="mockMrna", name_shuffle="nucleotides_method1", name_of_file="darnell_human_ViennaDuplex_features_negative")
    #
    # NegativeMockMrna.main(name_shuffle='denucleotides')
    # full_pipline(name_of_method="mockMrna", name_shuffle="denucleotides_method1", name_of_file="darnell_human_ViennaDuplex_features_negative")
    #

    ####### mockMrna 2 approach ######
    #5266684
    # from pipeline_steps.duplex_step import duplex as duplex_positive
    # from generate_interactions.mockMrna import NegativeMockMrna2
    # NegativeMockMrna2.main(name_shuffle='nucleotides')
    # full_pipline(name_of_method="mockMrna", name_shuffle="nucleotides_method2", name_of_file="darnell_human_ViennaDuplex_features_negative", duplex=duplex_positive)
    # NegativeMockMrna2.main(name_shuffle='denucleotides')
    # full_pipline(name_of_method="mockMrna", name_shuffle="denucleotides_method2", name_of_file="darnell_human_ViennaDuplex_features_negative", duplex=duplex_positive)

    ####### mockMrna 3 approach ######
    #5266685
    from pipeline_steps.duplex_step import duplex as duplex_positive
    from generate_interactions.mockMrna import NegativeMockMrna3
    NegativeMockMrna3.main(name_shuffle='nucleotides')
    full_pipline(name_of_method="mockMrna", name_shuffle="nucleotides_method3", name_of_file="darnell_human_ViennaDuplex_features_negative", duplex=duplex_positive)
    NegativeMockMrna3.main(name_shuffle='denucleotides')
    full_pipline(name_of_method="mockMrna", name_shuffle="denucleotides_method3", name_of_file="darnell_human_ViennaDuplex_features_negative", duplex=duplex_positive)
