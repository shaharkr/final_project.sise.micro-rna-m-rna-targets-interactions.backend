from dal.db_connection import cache, db
from dal.db_connection import Interaction, DataSet
from flask import Response, stream_with_context
import csv
from configurator import Configurator
import uuid
from sqlalchemy import or_


@cache.memoize(timeout=12000)
def get_interactions_for_download(data_sets_ids, seed_families, mirna_ids,
                     mirna_seqs, sites, gene_ids, regions, cur_file_num):
    interactions = None
    try:
        filters = [Interaction.data_set_id.in_(data_sets_ids) if data_sets_ids else True,
           Interaction.seed_family.in_(seed_families) if seed_families else True,
           Interaction.mirna_id.in_(mirna_ids) if mirna_ids else True,
           Interaction.mirna_sequence.in_(mirna_seqs) if mirna_seqs else True,
           Interaction.site.in_(sites) if sites else True,
           Interaction.Gene_ID.in_(gene_ids) if gene_ids else True,
           Interaction.region.in_(regions) if regions else True]
        interactions = Interaction.query.filter(*filters).all()
        file_path = create_csv_from_search(interactions, cur_file_num)
    except Exception as e:
        print(f'dal failed to get interactions. error: {str(e)}')
    return file_path

@cache.memoize(timeout=12000)
def get_interactions_gneral_search_for_download(query_string, file_name):
    interactions = []
    if query_string is None or query_string == '':
        return interactions
    try:
        interactions = db.session.query(Interaction).filter(or_(
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
        )).all()
        file_path = create_csv_from_search(interactions, file_name)
    except Exception as e:
        print(f'dal failed to get general interactions. error: {str(e)}')
    return file_path

def create_csv_from_search(interactions, file_name):
    # Check if the interactions list is not empty
    if interactions:
        try:
            # Set the file name and location
            path_prefix = Configurator().get_path_prefix_to_save_new_csv()
            file_path = f"{path_prefix}\\temp_{file_name}.csv"
            # Open the file in write mode and create a CSV writer object
            with open(file_path, mode='w', newline='') as csv_file:
                writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                # Write the header row
                writer.writerow(['index', 'datasetId', 'miRnaSeq', 'seedFamily', 'site', 'region', 'start', 'end', 'mrnaBulge', 
                                'mrnaInter', 'mirInter', 'mirBulge', 'energyMefDuplex', 'geneId'])
                # Loop through the interactions and write each row to the CSV file
                for interaction in interactions:
                    writer.writerow([interaction.index, interaction.data_set_id, interaction.mirna_sequence, interaction.seed_family, 
                                    interaction.site, interaction.region, interaction.start, interaction.end, interaction.mrna_bulge, 
                                    interaction.mrna_inter, interaction.mir_inter, interaction.mir_bulge, interaction.Energy_MEF_Duplex, 
                                    interaction.Gene_ID])
            return file_path
        except Exception as e:
            print(f'dal failed to get interactions. error: {str(e)}')
        
    else:
        print('Interactions list is empty.')

# if path=None - download a whole dataset
def download_data(data_set_id, path):
    # set the file path and content type
    file_path = path
    file_name = "search_data"
    if not path:
        path_prefix = Configurator().get_path_prefix_of_dataset_location()
        data_set = DataSet.query.filter_by(id=data_set_id).all()
        file_name = data_set[0].name
        file_path = f"{path_prefix}\\{file_name}.csv"
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
        response.headers['Content-Disposition'] = f'attachment; filename={file_name}.csv'
        return response
    except Exception as e:
        print(f'app failed to get general interactions. error: {str(e)}')
        return None
    

def get_unique_file_name():
    return str(uuid.uuid4())

