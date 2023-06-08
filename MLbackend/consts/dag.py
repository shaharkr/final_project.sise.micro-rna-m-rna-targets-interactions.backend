from airflow.utils.dates import days_ago
from datetime import timedelta

FILE_NAMES = ['dynamic_fly',
              'cattle',
              'human_mapping',
              'darnell_mouse',
              'darnell_human',
              'unambiguous_human',
              'unambiguous_celegans',
              'unambiguous_mouse',
              'pairing_beyond',
              'qclash_melanoma_human']

SPLUNK_DICT = {'dynamic_fly':("dynamic_fly", "fly"),
              'cattle':('cattle', 'cattle'),
              'human_mapping':('human_mapping', 'human'),
              'darnell_mouse': ('darnell', 'mouse'),
              'darnell_human': ('darnell', 'human'),
              'unambiguous_human': ('unambiguous', 'human'),
              'unambiguous_celegans': ('unambiguous', 'celegans'),
              'unambiguous_mouse': ('unambiguous', 'mouse'),
              'pairing_beyond':('pairing_beyond', 'celegans'),
              'qclash_melanoma_human':('qclash_melanoma', 'human')
               }


# Dag Names
MIRBASE_DOWNALOD = "mirbase_downalod"
FULL_PIPELINE_NAME = "full_pipeline"
DYNAMIC_FLY_GAMBIAE = "Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae"
HUMAN_MAPPING = "Mapping_the_Human_miRNA_Interactome_by_CLASH_Reveals_Frequent_Noncanonical_Binding"
CATTLE = "Global_mapping_of_miRNA-target_interactions_in_cattle_Bos_taurus"
DARNELL = "Darnell_miRNA_targe_chimeras_reveal"
UNAMBIGUOUS = "Unambiguous_Identification_of_miRNA_Target_Site_Interactions"
PAIRING_BEYOND = "Pairing_beyond_the_Seed_Supports_MicroRNA_Targeting_Specificity"
QCLASH_MELANOMA = "Cross-Linking_Ligation_and_Sequencing_of_Hybrids_qCLASH_Reveals_an_Unpredicted_miRNA_Targetome_in_Melanoma_Cells"

default_args = dict(owner='Gilad Ben Or', depends_on_past=False, start_date=days_ago(31),
                    email=['benorgi@post.bgu.ac.il'], email_on_failure=False, email_on_retry=False, retries=1,
                    retry_delay=timedelta(minutes=1))

#schedule_interval='None')
