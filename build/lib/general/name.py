import gffutils
import pandas as pd
import numpy as np


def gene_id_to_name(db):
    '''
    shamelessly taken from gscripts.region_helpers
    '''
    db = gffutils.FeatureDB(db)
    genes = db.features_of_type('gene')
    gene_name_dict = {}
    for gene in genes:
        gene_id = gene.attributes['gene_id'][0] if type(gene.attributes['gene_id']) == list else gene.attributes[
            'gene_id']
        try:
            gene_name_dict[gene_id] = gene.attributes['gene_name'][0]
        except KeyError:
            print("Warning. Key not found for {}".format(gene))
    return gene_name_dict


def gene_id_to_type(db):
    db = gffutils.FeatureDB(db)
    genes = db.features_of_type('gene')
    gene_type_dict = {}
    for gene in genes:
        gene_id = gene.attributes['gene_id'][0] if type(gene.attributes['gene_id']) == list else gene.attributes[
            'gene_id']
        try:
            gene_type_dict[gene_id] = gene.attributes['gene_type'][0]
        except KeyError:
            pass
    return gene_type_dict


def get_annotation_df(db, namecol='Name', genecol='Geneid'):
    df = pd.DataFrame.from_dict(gene_id_to_name(db), orient='index')
    df.columns = [namecol]
    df.index.rename(genecol, inplace=True)
    return df


def map_name_to_df(db, df, as_index=False):
    mapping = get_annotation_df(db)
    merged = pd.merge(mapping, df, how='right', left_index=True, right_index=True)
    if as_index:
        merged.set_index(mapping.columns[0], inplace=True)
    return merged
