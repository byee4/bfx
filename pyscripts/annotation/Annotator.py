import pybedtools
from tqdm import trange
import gff_helpers as g
import region_helpers as r
import gffutils
from collections import defaultdict
import os

PRIORITY = ['CDS', 'UTR', 'GENE', 'EXON', 'TRANSCRIPT', 'START_CODON', 'STOP_CODON']

HASH_VAL = 1000000
MAXVAL = 1000000000
MINVAL = 0

class Annotator():
    """

    class to get genomic features from gffutils _db

    """

    def __init__(self, db_file, chroms=[]):
        """
        
        :param db_file: 
        :param chroms: 
        """

        self._db = gffutils.FeatureDB(db_file)

        self._featuretypes = self._db.featuretypes()
        self._geneid_to_name_dict = self._gene_id_to_name()
        if len(chroms) != 0:
            self._chromosomes = set(chroms)
        else:
            self._chromosomes = self._chromosome_set()

        self.features_dict = self._hash_features()
        self.cds_dict = self._get_all_cds_dict()
        # self.genes_dict = g.gene_id_to_transcript

    def _chromosome_set(self):
        """
        Returns the set of chromosomes that exist in a database.
        
        :return: 
        """
        ret = self._db.execute("SELECT seqid FROM features").fetchall()
        all_chromosomes = [r['seqid'] for r in ret]
        return set(all_chromosomes)

    def _hash_features(self):
        """
        hashes features by position.
        :return features_dict : collections.defaultdict()
            dictionary of features{[chrom, pos/HASH_VAL, strand] : feature_list}
        """

        features_dict = defaultdict(list)
        progress = trange(len(self._chromosomes), leave=False, desc='Build Location Index')
        for chrom in self._chromosomes:
            for element in self._db.region(seqid=chrom):
                start = int(element.start / HASH_VAL)
                end = int(element.end / HASH_VAL)
                for i in range(start, end+1):
                    features_dict[chrom, i, element.strand].append(element)
            progress.update(1)
        return features_dict

    def _get_all_cds_dict(self):
        """
        For every cds-annotated transcript id (ENST), return a 
        dictionary containing the cds start and end vals for that transcript.
        
        :return cds_dict : defaultdict{transcript:{'start':START, 'end':END}} 
        """
        cds_dict = defaultdict(lambda: {'start': MAXVAL, 'end': MINVAL})
        for cds_feature in self._db.features_of_type('CDS'):
            for transcript_id in cds_feature.attributes['transcript_id']:
                if cds_feature.start <= cds_dict[transcript_id]['start']:
                    cds_dict[transcript_id]['start'] = cds_feature.start
                if cds_feature.end >= cds_dict[transcript_id]['end']:
                    cds_dict[transcript_id]['end'] = cds_feature.end
        return cds_dict

    def _classify_utr(self, utr_feature):
        """
        Given a list of features, return a dictionary of 
        gene_id : {start: cds_start, end: cds_end} sites

        :param features : list[gffutils.Feature]
        :return cds_ends : dict
        """
        three_prime_utr = False
        five_prime_utr = False

        for transcript_id in utr_feature.attributes['transcript_id']:
            if utr_feature.strand == '+':
                if self.cds_dict[transcript_id]['start'] > utr_feature.end:
                    five_prime_utr = True
                if self.cds_dict[transcript_id]['end'] < utr_feature.start:
                    three_prime_utr = True
            elif utr_feature.strand == '-':
                if self.cds_dict[transcript_id]['start'] > utr_feature.end:
                    three_prime_utr = True
                if self.cds_dict[transcript_id]['end'] < utr_feature.start:
                    five_prime_utr = True

        if five_prime_utr and three_prime_utr:
            return 'three_and_five_prime_utr'
        elif five_prime_utr:
            return 'five_prime_utr'
        elif three_prime_utr:
            return 'three_prime_utr'
        else:
            return 'unclassified_utr'

    def _gene_id_to_name(self):
        """
        Returns a dictionary containing a gene_id:name translation
        Note: may be different if the 'gene_id' or 'gene_name' 
        keys are not in the source GTF file
        (taken from gscripts.region_helpers)
        
        :return gene_name_dict : dict
            dict of {gene_id : gene_name}
        """

        genes = self._db.features_of_type('gene')
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

    def get_all_overlapping_features_from_query(self, chrom, qstart, qend, strand):
        """
        Given a query location (chr, start, end), return all features that 
        overlap by at least one base. Functions similarly to gffutils db.region(),
        but uses the pre-hashed self.features_dict to greatly speed things up.

        :param chrom : string
        :param qstart : int
        :param qend : int
        :param strand : string
        :return features : list
            list of gffutils.Feature objects.
        """
        features = []
        start_key = int(qstart / HASH_VAL)
        end_key = int(qend / HASH_VAL)
        for i in range(start_key, end_key + 1):
            for feature in self.features_dict[chrom, i, strand]:
                if qstart <= feature.start and qend >= feature.end:  # feature completely contains query
                    features.append(feature)
                elif qstart >= feature.start and qend <= feature.end:  # query completely contains feature
                    features.append(feature)
                elif qstart <= feature.start and qend >= feature.start:  # feature partially overlaps (qstart < fstart < qend)
                    features.append(feature)
                elif qstart <= feature.end and qend >= feature.end:  # feature partially overlaps (qstart < fend < qend)
                    features.append(feature)
        return features




    def prioritize(self, features, v=True, priority=PRIORITY):
        """
        Groups each transcript feature into genes.
        """
        hi = []  # contains highest priority genic feature for each gene
        genes = defaultdict(list)
        for feature in features:
            for gene_id in feature.attributes['gene_id']:
                genes[gene_id].append(feature)
        for gene, features in genes.iteritems():
            features.sort(key=lambda x: priority.index(x.featuretype.upper()))
        for gene_id, overlapping_features in genes.iteritems():
            if v == True:
                print(gene_id)
                for feature in overlapping_features:
                    print(feature.featuretype, feature.attributes['transcript_id'], feature.attributes['gene_id'], feature.start, feature.end)
            first_priority = overlapping_features[0]
            if first_priority.featuretype == 'UTR':

                annotation = self._classify_utr(first_priority)
            else:
                annotation = first_priority.featuretype.replace('gene','INTRON')
            hi.append('{}|{}'.format(annotation, gene_id))
        return hi