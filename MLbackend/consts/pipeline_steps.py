from consts.global_consts import ROOT_PATH
from numpy import int64

READ_PATH = ROOT_PATH / "data/pipeline_steps/read"
MIRNA_SEQ_PATH = ROOT_PATH / "data/pipeline_steps/mirna_sequence"
SITE_PATH = ROOT_PATH / "data/pipeline_steps/site"
REGION_PATH = ROOT_PATH / "data/pipeline_steps/region"
CONCAT_BLAST = ROOT_PATH / "data/pipeline_steps/concat_blast"
NORMALIZATION_PATH = ROOT_PATH / "data/pipeline_steps/normalization_final"


# NORMALIZATION_COLUMNS = ['key', 'paper name', 'organism', 'miRNA ID', 'miRNA sequence', 'seed_family',
#                          'site',
#                          'region', 'paper region',
#                          'start', 'end', 'sequence', 'Gene_ID',
#                          'identity', 'coverage',
#                          'valid_row']

NORMALIZATION_COLUMNS = ['key', 'paper name', 'organism', 'miRNA ID', 'miRNA sequence', 'seed_family',
                         'site',
                         'region',
                         'start', 'end', 'full_mrna', 'Gene_ID', 'valid_row']




