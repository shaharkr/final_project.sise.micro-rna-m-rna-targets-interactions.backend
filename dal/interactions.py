from dal.db_connection import db, cache
from dal.db_connection import Interaction, DataSet
from sqlalchemy import or_, and_
from sqlalchemy import text
from flask import Response, stream_with_context
import io
import csv
from configurator import Configurator


@cache.memoize(timeout=12000)
def get_interactions(data_sets_ids, seed_families, mirna_ids,
                     mirna_seqs, site_types, gene_ids, regions):
    interactions = []
    try:
        filters = [Interaction.data_set_id.in_(data_sets_ids) if data_sets_ids else True,
                  Interaction.seed_family.in_(seed_families) if seed_families else True,
                  Interaction.mirna_id.in_(mirna_ids) if mirna_ids else True,
                  Interaction.mirna_sequence.in_(mirna_seqs) if mirna_seqs else True,
                  Interaction.Gene_ID.in_(gene_ids) if gene_ids else True,
                  Interaction.region.in_(regions) if regions else True]
        site_types_filters = None
        if site_types:
            site_types_filters = []
            if 'canonical' in site_types:
                site_types_filters.append(Interaction.Seed_match_canonical == True)
            if 'noncanonical' in site_types:
                site_types_filters.append(Interaction.Seed_match_noncanonical == True)
            if 'other' in site_types:
                site_types_filters.append(and_(
                    Interaction.Seed_match_canonical == False,
                    Interaction.Seed_match_noncanonical == False
                ))
            filters.append(or_(*site_types_filters))
        results = Interaction.query.filter(*filters).limit(750).all()
        interactions = create_interactions_list(results)
    except Exception as e:
        print(f'dal failed to get interactions. error: {str(e)}')
    return interactions


@cache.memoize(timeout=12000)
def get_interactions_general_search(query_string):
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
            interactions.append({"interactionId": interaction.id,
                                "index": interaction.index,
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


@cache.memoize(timeout=12000)
def get_interaction_id_data(interaction_id):
    interaction_id_data = {}
    try:
        results = Interaction.query.filter_by(id=interaction_id)
        interaction_id_data = create_interaction_outer_data_object(results[0])
    except Exception as e:
        print(f'dal failed to get general interactions. error: {str(e)}')
    return interaction_id_data


def create_interaction_outer_data_object(interaction_id_data):
    mi_rna_data_dict = {
            "miRnaId": interaction_id_data.mirna_id,
            "miRnaSequence": interaction_id_data.mirna_sequence,
            "seedFamily": interaction_id_data.seed_family
        }
    sequence_url = create_sequence_url(interaction_id_data.Gene_ID, interaction_id_data.organism)
    m_rna_data_dict = {
            "region": interaction_id_data.region,
            "geneId": interaction_id_data.Gene_ID,
            "geneName": "temp name",  # TODO: need to add this column to DB
            "sequenceUrl": sequence_url,
            "startSite": interaction_id_data.start,
            "endSite": interaction_id_data.end
        }
    
    duplex_structure = create_duplex_structure(interaction_id_data.mrna_bulge,
                                               interaction_id_data.mrna_inter,
                                               interaction_id_data.mir_inter,
                                               interaction_id_data.mir_bulge)
    interaction_inner_data_dict = {
            "interactionId": interaction_id_data.id,
            "organismName": interaction_id_data.organism,
            "dataSource": interaction_id_data.paper_name,
            "duplexStructure": duplex_structure,
            "energyMefDuplex": interaction_id_data.Energy_MEF_Duplex,
            "mRnaDistToEnd": interaction_id_data.MRNA_Dist_to_end,
            "mRnaDistToStart": interaction_id_data.MRNA_Dist_to_start,
            "seedMatchCanonical": interaction_id_data.Seed_match_canonical,
            "seedMatchNonCanonical": interaction_id_data.Seed_match_noncanonical,
            "seedMatchStart": interaction_id_data.Seed_match_start
        }
    
    interaction_outer_data_object = {
        "miRnaData": mi_rna_data_dict,
        "mRnaData": m_rna_data_dict,
        "interactionInnerData": interaction_inner_data_dict
        }
    return interaction_outer_data_object


def create_duplex_structure(mrna_bulge, mrna_inter, mir_inter, mir_bulge):
    mrna_bulge_prefix = "mRNA  3’ "
    mrna_bulge_sufix = " 5’"
    mrna_inter_prefix = "         "
    mrna_inter_sufix = "   "
    mir_inter_prefix = "         "
    mir_inter_sufix = "   "
    mir_bulge_prefix = "miRNA 5’ "
    mir_bulge_sufix = " 3’"
    return f'{mrna_bulge_prefix}{mrna_bulge}{mrna_bulge_sufix}\n{mrna_inter_prefix}{mrna_inter}{mrna_inter_sufix}\n{mir_inter_prefix}{mir_inter}{mir_inter_sufix}\n{mir_bulge_prefix}{mir_bulge}{mir_bulge_sufix}'


def create_sequence_url(gene_id, organism):
    confg = Configurator()
    relevant_gene_id = gene_id.split('|')[0]
    ensamble_org_name = confg.get_ensamble_organisms_names_dict()[organism]
    url_string = confg.get_ensambl_url().format(organism_name=ensamble_org_name, gene_id=relevant_gene_id)
    return url_string

  
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
    path_prefix = Configurator().get_path_prefix_of_dataset_location()
    data_set = DataSet.query.filter_by(id=data_set_id).all()
    file_name = data_set[0].name
    if not path:
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
