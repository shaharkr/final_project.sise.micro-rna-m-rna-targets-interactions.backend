from dal.db_connection import db, cache
from dal.db_connection import Interaction
# from main import app
from sqlalchemy import or_
from sqlalchemy import text
from flask import  make_response
import io
import csv


# def init_base(app):
#     Base = automap_base()
#     with app.app_context():
#         Base.prepare(autoload_with=db.engine, reflect=True)
#     return Base

@cache.memoize(timeout=12000)
def get_interactions(data_sets_ids, seed_families, mirna_ids,
                     mirna_seqs, sites, gene_ids, regions):
    interactions = []
    try:
        filters = [Interaction.data_set_id.in_(data_sets_ids) if data_sets_ids else True,
           Interaction.seed_family.in_(seed_families) if seed_families else True,
           Interaction.mirna_id.in_(mirna_ids) if mirna_ids else True,
           Interaction.mirna_sequence.in_(mirna_seqs) if mirna_seqs else True,
           Interaction.site.in_(sites) if sites else True,
           Interaction.Gene_ID.in_(gene_ids) if gene_ids else True,
           Interaction.region.in_(regions) if regions else True]
        results = Interaction.query.filter(*filters).limit(750).all()
        interactions = create_interactions_list(results)
    except Exception as e:
        print(f'dal failed to get interactions. error: {str(e)}')
    return interactions


@cache.memoize(timeout=12000)
def get_interactions_gneral_search(query_string):
    interactions = []
    if query_string is None or query_string == '':
        return interactions
    try:
        results = db.session.query(Interaction).filter(or_(
            Interaction.mirna_id.like(f'%{query_string}%'),
            Interaction.mirna_sequence.like(f'%{query_string}%'),
            Interaction.seed_family.like(f'%{query_string}%'),
            Interaction.site.like(f'%{query_string}%'),
            Interaction.region.like(f'%{query_string}%'),
            Interaction.mrna_bulge.like(f'%{query_string}%'),
            Interaction.mrna_inter.like(f'%{query_string}%'),
            Interaction.mir_inter.like(f'%{query_string}%'),
            Interaction.mir_bulge.like(f'%{query_string}%'),
            Interaction.Gene_ID.like(f'%{query_string}%')
        )).limit(750).all()
        interactions = create_interactions_list(results)
    except Exception as e:
        print(f'dal failed to get general interactions. error: {str(e)}')
    return interactions


def create_interactions_list(results):
    interactions = []
    for interaction in results:
            interactions.append({"index": interaction.index,
                                "datasetId": interaction.data_set_id,
                                "miRnaId": interaction.mirna_id,
                                "miRnaSeq": interaction.mirna_sequence,
                                "seedFamily": interaction.seed_family,
                                "site": interaction.site,
                                "region": interaction.region,
                                "start": interaction.start,
                                "end": interaction.end,
                                "mrnaBulge": interaction.mrna_bulge,
                                "mrnaInter": interaction.mrna_inter,
                                "mirInter": interaction.mir_inter,
                                "mirBulge": interaction.mir_bulge,
                                "energyMefDuplex": interaction.Energy_MEF_Duplex,
                                "geneId": interaction.Gene_ID
                            })
    return interactions


def downlod_dataset(data_set_id):
    # Base = init_base(app)
    # mirna_mrna_interactions = Base.classes.mirna_mrna_interactions
    # for class_ in Base.classes:
    #     print(class_)
    # interactions_dict = []
    try:
        page_size = 10000
        with db.engine.connect() as conn:
            result = conn.execute(text("SELECT * FROM mirna_mrna_interactions WHERE data_set_id ="+ str(data_set_id)))
        # Paginate the result set
        offset = 0
        rows = result.fetchmany(page_size)

        # Write the query result to a CSV file in memory
        csv_data = io.StringIO()
        writer = csv.writer(csv_data)
        writer.writerow([column[0] for column in result.cursor.description])
        while rows:
            for row in rows:
                writer.writerow(row)

            # Fetch the next page
            offset += page_size
            rows = result.fetchmany(page_size)

        # Create the response object
        csv_data.seek(0)
        response = make_response(csv_data.getvalue().encode('utf-8'))
        response.headers['Content-Disposition'] = 'attachment; filename=data.csv'
        response.headers['Content-Type'] = 'text/csv'
        response.headers['Content-Encoding'] = 'utf-8'

        return response
        # res = db.session.query(mirna_mrna_interactions).filter_by(data_set_id=data_set_id).all()
        # for inter in res:
        #     interactions_dict.append(inter)

        # interactions_dict = [interaction.to_dict() for interaction in results]
    except Exception as e:
        print(f'dal failed to get general interactions. error: {str(e)}')

    