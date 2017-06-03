import pandas as pd
import gffutils
import os
from collections import defaultdict

GTF_NAMES = ['chrom','source','feature_type','start','end','.','strand','.','attributes']


def get_feature_type_set(gtf_file):
    """
    from a GTF file, extract the set of feature_types
    (feature_types is the third column, normally)
    This might be useful for figuring out the priority for annotation.
    
    Parameters
    ----------
    gtf_file

    Returns
    -------

    """
    gtf_df = pd.read_table(
        gtf_file,
        names=GTF_NAMES,
        comment='#'
    )
    return set(gtf_df['feature_type'])


def get_attribute_type_set(gtf_file, attribute_type):
    """
    from a GTF file, extract the set of attribute_types
    (attribute_types is one of those fields contained within the 9th column)
    This might be useful for figuring out the priority for annotation.
    
    Parameters
    ----------
    gtf_file : basestring
    attribute_type : basestring

    Returns
    -------

    """

    gtf_df = pd.read_table(
        gtf_file,
        names=GTF_NAMES,
        comment='#'
    )
    regex_filter = '{} \"([\w\s\d -]+)\"'.format(attribute_type)
    return set(gtf_df['attributes'].str.extract(regex_filter, expand=False))


def build_db(annotation_file, annotation_db_file, force=False):
    """
    builds a database file from a GTF/GFF annotation file
    
    Parameters
    ----------
    annotation_file : basestring
    annotation_db_file : basestring
    force: bool
        raises exception if db (annotation_db_file) exists
    Returns
    -------
    0 if successful, 1 if not.
    """
    try:
        gffutils.create_db(
            annotation_file, dbfn=annotation_db_file, force=force,
            keep_order=True, merge_strategy='merge',
            sort_attribute_values=True
        )
        return 0
    except Exception as e:
        print(e)
        return 1


def gene_id_to_name(db):
    """
    Returns a dictionary containing a gene_id:name translation
    Note: may be different if the 'gene_id' or 'gene_name' 
    keys are not in the source GTF file
    (taken from gscripts.region_helpers)
    
    Parameters
    ----------
    db : database

    Returns
    -------

    """

    genes = db.features_of_type('gene')
    gene_name_dict = {}
    for gene in genes:
        gene_id = gene.attributes['gene_id'][0] if type(gene.attributes['gene_id']) == list else gene.attributes[
            'gene_id']
        try:
            gene_name_dict[gene_id] = gene.attributes['gene_name'][0]
        except KeyError:
            print(gene.attributes.keys())
            print("Warning. Key not found for {}".format(gene))
            return 1
    return gene_name_dict


def gene_name_to_id(db):
    """
    given a gene name, returns a list of associated Gene IDs (one-to-many)
    
    Parameters
    ----------
    db : sqlite3 database

    Returns
    -------
    gene_name_dict : dict{list}
        dictionary of gene_name : [gene_ids]
    """

    genes = db.features_of_type('gene')
    gene_name_dict = defaultdict(list)
    for gene in genes:
        try:
            gene_name_dict[gene.attributes['gene_name'][0]].append(gene.attributes['gene_id'][0])
        except KeyError as e:
            print("Warning. Key not found for {}".format(gene))
            return 1
    return gene_name_dict


def gene_name_to_transcript(db):
    """
    given a gene name, returns a list of associated transcript IDs (one-to-many)
    
    Parameters
    ----------
    db : sqlite3 database

    Returns
    -------
    gene_name_dict : dict{list}
        dictionary of gene_name : [transcript_ids]
    """

    genes = db.features_of_type('transcript')
    gene_name_dict = defaultdict(list)
    for gene in genes:
        try:
            gene_name_dict[gene.attributes['gene_name'][0]].append(gene.attributes['transcript_id'][0])
        except KeyError as e:
            print("Warning. Key not found for {}".format(gene))
            return 1
    return gene_name_dict


def id_to_exons(db, identifier):
    """
    takes the gene or transcript id and returns exon positions
    
    Parameters
    ----------
    db : sqlite3 database
    identifier : basestring
        gene_id or gene_name

    Returns
    -------
    exons : list
        list of exons<Feature> corresponding to the gene
    """

    exons = []
    for i in db.children(identifier, featuretype='exon', order_by='start'):
        exons.append(i)
    return exons


def position_to_features(db, chrom, start, end, strand='', completely_within=True):
    """
    
    takes a coordinate and returns all the features overlapping 
    (either completely contained or partially overlapping the region).
    
    Parameters
    ----------
    db : sqlite3 database
    chrom : basestring
    start : int
    end : int
    strand : basestring
    completely_within : bool
        True if the features returned must be completely contained
        within the region. False if the features need only to be
        partially overlapping the region.
        
    Returns
    -------
    region_list: list
        list of Features corresponding to overlapping/contained
        features intersecting a region.
    """

    if strand == '+' or strand == '-':
        return list(
            db.region(
                region=(chrom, start, end), strand=strand, completely_within=completely_within
            )
        )
    else:
        return list(
            db.region(
                region=(chrom, start, end), completely_within=completely_within
            )
        )


def bedtool_to_features(db, interval, completely_within):
    """
    
    takes a coordinate and returns all the features overlapping 
    (either completely contained or partially overlapping the region).
    
    Parameters
    ----------
    db : sqlite3 database
    interval : pybedtools.Interval
        interval object
    completely_within : bool
        True if the features returned must be completely contained
        within the region. False if the features need only to be
        partially overlapping the region.
        
    Returns
    -------
    region_list: list
        list of Features corresponding to overlapping/contained
        features intersecting a region.
    """
    return position_to_features(
        db,
        interval.chrom,
        interval.start,
        interval.end,
        interval.strand,
        completely_within
    )

