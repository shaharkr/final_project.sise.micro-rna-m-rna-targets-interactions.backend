from dal.db_connection import GeneralStats, Interaction, DataSet, Organism
from dal.db_connection import db, cache
from sqlalchemy import func, distinct
import pandas as pd


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


@cache.memoize(timeout=12000)
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

