MAXVAL=10000000000
MINVAL=0

from collections import defaultdict

def get_cds_starts(features):
    """
    Given a list of features, return a dictionary of 
    transcript_id : cds_start sites
    
    :param features : list[gffutils.Feature]
    :return cds_ends : dict
    """
    cds = [f for f in features if f.featuretype == 'CDS']

    cds_starts = defaultdict(lambda: MAXVAL)

    for c in cds:
        for transcript in c.attributes['transcript_id']:
            if c.start <= cds_starts[transcript]:
                cds_starts[transcript] = c.start


def get_cds_ends(features):
    """
    Given a list of features, return a dictionary of 
    transcript_id : cds_end sites

    :param features : list[gffutils.Feature]
    :return cds_ends : dict
    """
    cds = [f for f in features if f.featuretype == 'CDS']

    cds_ends = defaultdict(lambda: MINVAL)

    for c in cds:
        for transcript in c.attributes['transcript_id']:
            if c.end >= cds_ends[transcript]:
                cds_ends[transcript] = c.end

def get_cds_dict(features):
    """
    Given a list of features, return a dictionary of 
    transcript_id : {start: cds_start, end: cds_end} sites

    :param features : list[gffutils.Feature]
    :return cds_ends : dict
    """
    cds = [f for f in features if f.featuretype == 'CDS']

    cds_dict = defaultdict(lambda:{'start':MAXVAL,'end':MINVAL})

    for c in cds:
        for transcript in c.attributes['transcript_id']:
            if c.start <= cds_dict[transcript]['start']:
                cds_dict[transcript]['start'] = c.start
            if c.end >= cds_dict[transcript]['end']:
                cds_dict[transcript]['end'] = c.end
    return cds_dict


def get_utrs_list(features, cds_dict):
    """
    Given a list of features, and a dictionary of regions defining
    cds start and end sites, return a dictionary of 
    transcript_id : {start: cds_start, end: cds_end} sites

    :param features : list[gffutils.Feature]
    :return five_prime_utrs : list
    :return three_prime_utrs : list
    """
    three_prime_utrs = []
    five_prime_utrs = []

    utr = [f for f in features if f.featuretype == 'UTR']
    for u in utr:
        for transcript_id in u.attributes['transcript_id']:
            if u.strand == '+':
                if cds_dict[transcript_id]['start'] > u.end:
                    five_prime_utrs.append(u)
                if cds_dict[transcript_id]['end'] < u.start:
                    three_prime_utrs.append(u)
            elif u.strand == '-':
                if cds_dict[transcript_id]['start'] > u.end:
                    three_prime_utrs.append(u)
                if cds_dict[transcript_id]['end'] < u.start:
                    five_prime_utrs.append(u)
    return five_prime_utrs, three_prime_utrs




def classify(features):
    cds_dict = get_cds_dict(features) # get all associated cds
    regions = defaultdict(list)

    regions['cds'] = [f for f in features if f.featuretype == 'CDS']
    regions['utr5'], regions['utr3'] = get_utrs_list(features, cds_dict)
    regions['utr53'] = list(set(regions['utr5']).intersection(set(regions['utr3'])))
    regions['exons'] = [f for f in features if f.featuretype == 'exon']
    return regions