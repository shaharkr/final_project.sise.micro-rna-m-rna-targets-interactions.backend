from dal.db_connection import db, cache
from dal.db_connection import Interaction
from sqlalchemy import or_
from sqlalchemy import text
from flask import Response, stream_with_context
import io
import csv


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


def download_search_data(data_set_id):
    page_size = 10240
    with db.engine.connect() as conn:
        result = conn.execute(text("SELECT * FROM mirna_mrna_interactions WHERE data_set_id =" + str(data_set_id)))
    
    def generate():
        offset = 0
        rows = result.fetchmany(page_size)
        csv_data = io.StringIO()
        writer = csv.writer(csv_data)
        writer.writerow([column[0] for column in result.cursor.description])
        
        while rows:
            for row in rows:
                writer.writerow(row)
                data = csv_data.getvalue()
                if len(data) > 10240: # if csv data exceeds 10KB
                    yield data
                    csv_data = io.StringIO() # reset the StringIO object
                    writer = csv.writer(csv_data)
                    writer.writerow([column[0] for column in result.cursor.description])
            offset += page_size
            rows = result.fetchmany(page_size)

        yield csv_data.getvalue() # yield any remaining csv data
        
    content_type = 'text/csv'
    return Response(stream_with_context(generate()), mimetype=content_type)


def download_data(data_set_id, path):
    # set the file path and content type
    file_path = path
    if not path:
        file_path = f"C:\\datasets\\{data_set_id}.csv"
    content_type = 'text/csv'

    # define a function that reads the file in 10KB chunks
    def generate():
        with open(file_path, 'rb') as f:
            while True:
                data = f.read(10240) # 10KB chunk size
                if not data:
                    break
                yield data
    
    try:
    # use the stream_with_context function to stream the response in chunks
        response = Response(stream_with_context(generate()), mimetype=content_type)
        response.headers['Content-Disposition'] = f'attachment; filename={data_set_id}.csv'
        return response
    except Exception as e:
        print(f'app failed to get general interactions. error: {str(e)}')
        return None

    