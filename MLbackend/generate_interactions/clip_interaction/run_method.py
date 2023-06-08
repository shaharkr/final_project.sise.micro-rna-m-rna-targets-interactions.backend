from generate_interactions.clip_interaction.mirna_files import run_mirna_generate_files
from generate_interactions.clip_interaction.mrna_files import run_mrna_generate_files
from generate_interactions.clip_interaction.generate import genetate_interactions_from_files

def run():
    # step 1 ---> generate miRNA files
    # run_mirna_generate_files()

    # step 2 ---> generate mRNA files

    run_mrna_generate_files()

    # step 3 ---> combine two files and get final file of interactions
    genetate_interactions_from_files()