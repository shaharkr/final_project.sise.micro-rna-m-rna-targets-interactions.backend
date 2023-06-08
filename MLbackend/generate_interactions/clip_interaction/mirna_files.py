from utils.utilsfile import read_csv, to_csv, get_wrapper
from consts.global_consts import ROOT_PATH, DATA_PATH,CLIP_PATH_DATA
import pandas as pd
from consts.mirna_utils import MIRBASE_FILE

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



def read_mirna_files():
    file_name = MIRBASE_FILE
    df_mirna = pd.read_csv(file_name)
    df_mirna_hsa = df_mirna[df_mirna['miRNA ID'].str.contains('hsa')]
    df_mirna_hsa['version'] = pd.to_numeric(df_mirna_hsa['version'])
    df_mirna_hsa = df_mirna_hsa.loc[df_mirna_hsa.groupby(['miRNA ID'])["version"].idxmax()]
    df_mirna_hsa['miRNA_ID_perfix'] = df_mirna_hsa['miRNA ID'].apply(lambda name: name.split("-")[2])
    return df_mirna_hsa


def mirna_pre_mirna_to_mature(pre_mirna, sequence_precosor, df_mirna_hsa):
    # print(pre_mirna)
    # if sequence_precosor == 'UCAGUUAUCACAGUGCUGAUG':
    #     print("f")
    pre_mirna_name = pre_mirna.split("-")[2]
    df_mirna_hsa['miRNA_ID_perfix'] = df_mirna_hsa['miRNA_ID_perfix'].apply(lambda s: s.lower())
    df_mirna_hsa_match = df_mirna_hsa[df_mirna_hsa['miRNA_ID_perfix'] == pre_mirna_name.lower()]

    df_mirna_hsa_match['max_substring_one_direct'] = df_mirna_hsa_match['miRNA sequence'].apply(lambda sequence_mature: substringFinder(sequence_mature, sequence_precosor))
    df_mirna_hsa_match['max_substring_second_direct'] = df_mirna_hsa_match['miRNA sequence'].apply(lambda sequence_mature: substringFinder(sequence_precosor, sequence_mature))


    df_mirna_hsa_match['overlap_len'] = df_mirna_hsa_match.apply(lambda x : max(len(x.max_substring_one_direct),len(x.max_substring_second_direct)), axis=1)
    max_val = max(df_mirna_hsa_match['overlap_len'])
    df_mirna_hsa_match = df_mirna_hsa_match[df_mirna_hsa_match['overlap_len'] == max_val]


    # check if the same mirna sequence were found
    df_mirna_hsa_match_last_version = df_mirna_hsa_match.loc[df_mirna_hsa_match.groupby(['miRNA sequence'])["version"].idxmax()]

    # value to return
    mature_mirna_ID = df_mirna_hsa_match_last_version['miRNA ID'].values[0]
    mature_mirna_sequence = df_mirna_hsa_match_last_version['miRNA sequence'].values[0]

    return mature_mirna_ID, mature_mirna_sequence



def creat_data_frame_list(read_list):

    name_cols = ["miRNA original", "miRNA ID", "sequence_original", "sequence"]
    df_list = []
    df_mirna_hsa = read_mirna_files()
    # read_list = read_list[:500]
    for ls in read_list:
        pre_mirna = ls[2]
        sequence = ls[9].replace('T', 'U')

        # don't save sequence that the len is greater then 25
        if len(sequence) > 25:
            continue

        mature_mirna_name, mature_mirna_sequence = mirna_pre_mirna_to_mature(pre_mirna, sequence, df_mirna_hsa)
        ls = [pre_mirna, mature_mirna_name, sequence,  mature_mirna_sequence]
        df_list.append(ls)
    df = pd.DataFrame(df_list, columns=name_cols, dtype=object)
    df = df.reset_index()
    print("size before group by:", df.shape)
    df_g = df.groupby(['miRNA ID'])
    df = df.loc[df_g["index"].idxmax()]
    df.reset_index(inplace=True, drop= True)
    print("size after group by:", df.shape)

    return df


def run_mirna_generate_files():
    clip_data_path = CLIP_PATH_DATA
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*result_mirna*"):
            # if "clip_1" not in str(file):
            #     continue
            read_list = []
            with open(file) as f:
                lines = f.readlines()
            for line in lines:
                if "@SQ" in line or "@HD" in line or "@PG" in line:
                    continue
                else:
                    split_line = line.split("\t")
                    read_list.append(split_line)

            # generate dataframe from the list
            df = creat_data_frame_list(read_list)
            name_file = "mirna.csv"
            path_df = clip_dir / name_file
            to_csv(df, path=path_df)


# run_mirna_generate_files()

# df = read_csv("/sise/vaksler-group/IsanaRNA/miRNA_target_rules/Isana/clip_10/mirna.csv")
# df_g = df.groupby(['miRNA ID'])
# print(df_g.ngroups)

