from utils.utilsfile import read_csv, to_csv, get_wrapper, get_subsequence_by_coordinates
from consts.global_consts import ROOT_PATH, DATA_PATH,CLIP_PATH_DATA
import pandas as pd
from consts.global_consts import BIOMART_PATH
from consts.mirna_utils import MIRBASE_FILE
from consts.global_consts import HUMAN_SITE_EXTENDED_LEN, SITE_EXTRA_CHARS


def read_mrna_files():
    # path = "/sise/vaksler-group/IsanaRNA/miRNA_target_rules/benorgi/pipeline/data/biomart/human_3utr.fa"
    # f = read_csv(path)
    file_name = BIOMART_PATH / "human_3utr.csv"
    df_human_3utr = pd.read_csv(file_name)
    df_human_3utr['Gene_ID'] = df_human_3utr['ID'].apply(lambda x: x.split("|")[0])

    # In cases where multiple UTRs exist per gene, we considered the
    # longest UTR.
    # df_human_3utr = df_human_3utr.loc[df_human_3utr.groupby(['Gene_ID'])["sequence length"].idxmax()]
    df_human_3utr = df_human_3utr.rename(columns={'sequence': 'full_mrna'})
    df_human_3utr['full_mrna'] = df_human_3utr['full_mrna'].apply(lambda seq: seq.replace('T', 'U'))
    return df_human_3utr

def complete_site_chars(start, end):
    len_site = end - start
    number_chars_add_one_side = 0
    if len_site < 75:
        number_chars_add = 75 - len_site
        number_chars_add_one_side = round(number_chars_add / 2)

    return number_chars_add_one_side

# print(complete_site_chars(50, 72))


def diff_letters(a,b):
    return sum(a[i] != b[i] for i in range(len(a)))

def find_full_mrna(df, df_mrna_hsa):
    # step 1 ---> find for each gen the full mrna by geneID+transcript
    df["key"] = df.reset_index().index
    df_new = pd.merge(df, df_mrna_hsa, on=['ID'])

    # step 2 ---> find for each read the longest transcript
    print("df before group:", df_new.shape)
    df_new = df_new.loc[df_new.groupby(['read_ID'])["sequence length"].idxmax()]
    print("df after group:", df_new.shape)

    # step 3 ---> find the new site on the full mrna
    df_new['len_site'] = df_new['site'].apply(lambda site: len(site))
    df_new['start'] = df_new['start'].astype(int)
    df_new['end'] = df_new['start'] + df_new['len_site'] - 1
    df_new["site_new"] = df_new.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "full_mrna", "start", "end",strand="+",
                                           extra_chars=0),axis=1)

    # step 4 ---> filer the mrna that doesn't match exactly to the original site
    print("data frame size before filter mismatch: ", df_new.shape)
    df_new["different"] = df_new.apply(func=get_wrapper(diff_letters, "site_new", "site"), axis=1)
    df_new = df_new[df_new["different"] == 0]
    print("data frame size before filter mismatch: ", df_new.shape)

    # step 4 ---> find the site with the extra chars
    # step 4a ---> find the number of chars that need in order to complete to 75 len

    df_new["number_chars_complete"] = df_new.apply(func=get_wrapper(complete_site_chars,
                                                       "start", "end"), axis=1)

    df_new['number_chars_complete'] = df_new['number_chars_complete'].astype(int)
    df_new['strand'] = "+"

    df_new["site_new"] = df_new.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                                       "full_mrna", "start", "end", "strand", "number_chars_complete"
                                                       ), axis=1)

    return df_new


def creat_data_frame_list(read_list):

    name_cols = ['read_ID',"ID", "site", 'start', 'mismatch']
    df_list = []
    for ls in read_list:
        read_ID = ls[0]
        mrna_ID = ls[2]
        site = ls[9].replace('T', 'U')
        strand = ls[1]
        start = ls[3]
        mismatch = ls[5]

        if strand != '0' and strand != '256':
            continue

        # don't save sequence that the len is greater then 25
        if len(site) < 40:
            continue

        # full_mrna, start, end = find_full_mrna(mrna_ID, site, df_mrna_hsa)
        item = [read_ID, mrna_ID, site, start, mismatch]
        df_list.append(item)
    df = pd.DataFrame(df_list, columns=name_cols, dtype=object)
    return df


def filter_interaction(df):
    df["key"] = df.reset_index().index
    print("dataframe before filter:", df.shape)
    df_new = df.loc[df.groupby(['read_ID'])["sequence length"].idxmax()]
    print("dataframe after filter:", df_new.shape)

    return df_new


def drop_duplicate(result):
    # Gene_ID ---> number Gen
    # ID ---> number Gen + number Transcript

    print("size before filter grop mrna:", result.shape)
    df_g = result.groupby(['ID', 'site_new'])
    result = result.loc[df_g["sequence length"].idxmax()]
    print("size after filter grop mrna:", result.shape)

    return result


def split_files_add_data():

    df_mrna_hsa = read_mrna_files()
    clip_data_path = CLIP_PATH_DATA
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*result_mrna*"):
            if "clip_3" not in str(file):
                continue
            read_list = []
            number_file = 0
            count = 0
            with open(file) as f:
                 for line in f:
                    if "@SQ" in line or "@HD" in line or "@PG" in line:
                        continue
                    else:
                        split_line = line.split("\t")
                        read_list.append(split_line)
                        count = count + 1

                    if count == 2000000:

                        # from list to dataframe
                        df = creat_data_frame_list(read_list)
                        df = find_full_mrna(df, df_mrna_hsa)
                        df = drop_duplicate(df)

                        name = str(number_file) + ".csv"
                        path_df = clip_dir / "split_file" / name
                        number_file = number_file + 1
                        to_csv(df, path_df)
                        count = 0
                        df._clear_item_cache()
                        read_list.clear()

                 # the buffer is close but there are more line
                 if count > 0:
                    # from list to dataframe
                    df = creat_data_frame_list(read_list)
                    df = find_full_mrna(df, df_mrna_hsa)
                    df = drop_duplicate(df)

                    name = str(number_file) + ".csv"
                    path_df = clip_dir / "split_file" / name
                    to_csv(df, path_df)
                    df._clear_item_cache()
                    read_list.clear()

def combine_files():

    clip_data_path = CLIP_PATH_DATA
    for clip_dir in clip_data_path.iterdir():
        clip_dir_new = clip_dir / "split_file"
        for f in clip_dir.glob("*result_mirna*"):
            if "clip_3" not in str(clip_dir):
                continue
            frames = []
            for file in clip_dir_new.glob("*csv*"):
                print("################################################")
                print(file.stem)
                df = read_csv(file)
                frames.append(df)

            result = pd.concat(frames)
            print(result.shape)
            result.reset_index(drop=True, inplace=True)
            result.reset_index(inplace=True)
            result = filter_interaction(result)
            result = drop_duplicate(result)
            name = "mrna" + ".csv"
            path_df = clip_dir / name
            to_csv(result, path_df)



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


def clean_mrna_files():
    clip_data_path = CLIP_PATH_DATA
    for clip_dir in clip_data_path.iterdir():
        index_remove = []

        for f in clip_dir.glob("*mrna.csv*"):
            mrna = read_csv(f)
            mrna.drop(columns=['index','key'], inplace=True)
            mrna.reset_index(drop=True, inplace=True)
            mrna.reset_index(inplace=True)
            df1_grouped = mrna.groupby(['ID'])

            print("***************************************************************************")
            for group_name, df_group in df1_grouped:
                for row_index1, row1 in df_group.iterrows():
                    if row_index1 in index_remove:
                        continue
                    start1 = row1['start']
                    end1 = row1['end']
                    site1 = row1['site_new']
                    for row_index2, row2 in df_group.iterrows():
                        start2 = row2['start']
                        end2 = row2['end']
                        site2 = row2['site_new']
                        if row_index1 != row_index2:
                            if row_index2 in index_remove:
                                continue
                            if start1 == start2 or start1 < start2:
                                if end1 == end2 or end2 < end1:
                                    index_remove.append(row_index2)
                                    continue

                            overlap_len = len(substringFinder(site1, site2))
                            site1_len = len(site1)
                            precent_cover = (overlap_len / site1_len) * 100
                            if precent_cover > 85:
                                index_remove.append(row_index2)
                                continue

            print(len(set(index_remove)))
            print("Size after filter:",mrna.shape[0] - len(set(index_remove)))
            # to = mrna.loc[mrna.index[index_remove]]
            new = mrna[(~mrna.index.isin(index_remove))]
            print(new.shape)

            path = f.parents[0] / "mrna_clean.csv"
            to_csv(new, path)



def run_mrna_generate_files():
    split_files_add_data()
    combine_files()
    clean_mrna_files()

# run_mrna_generate_files()



