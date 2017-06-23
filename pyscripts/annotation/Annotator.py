import pybedtools
from tqdm import trange
import gffutils
from collections import defaultdict
import copy
TRANSCRIPT_PRIORITY = ['CDS', '5UTR',  '3UTR', 'THREE_AND_FIVE_PRIME_UTR', 'EXON', 'GENE', 'NONCODING_EXON', 'NONCODING_INTRON', 'TRANSCRIPT', 'START_CODON', 'STOP_CODON']
GENE_PRIORITY = ['CDS', '5UTR', '3UTR', 'THREE_AND_FIVE_PRIME_UTR', 'GENE', 'EXON', 'NONCODING_EXON', 'NONCODING_INTRON', 'TRANSCRIPT', 'START_CODON', 'STOP_CODON']

# PRIORITY = ['CDS', '3UTR', '5UTR', 'THREE_AND_FIVE_PRIME_UTR', 'GENE', 'EXON', 'NONCODING_EXON', 'NONCODING_INTRON', 'TRANSCRIPT', 'START_CODON', 'STOP_CODON']


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
        progress = trange(8, desc='Initializing/creating defs')

        self._db = gffutils.FeatureDB(db_file)
        progress.update(1)
        self._featuretypes = self._db.featuretypes()
        progress.update(1)
        self._geneid_to_name_dict = self._gene_id_to_name()
        progress.update(1)
        if len(chroms) != 0:  # if use specific chromosomes, otherwise hash all
            self._chromosomes = set(chroms)
        else:
            self._chromosomes = self._chromosome_set()

        self.features_dict = self._hash_features()
        progress.update(1)
        self.exons_dict = self._get_all_exons_dict()
        progress.update(1)
        self.transcripts_dict = self._get_all_transcripts_dict()
        progress.update(1)
        self._update_introns()
        progress.update(1)
        self.cds_dict = self._get_all_cds_dict()
        progress.update(1)
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
        progress = trange(
            len(self._chromosomes),
            leave=False,
            desc='Build Location Index'
        )
        for chrom in self._chromosomes:
            for element in self._db.region(seqid=chrom):
                start = int(element.start / HASH_VAL)
                end = int(element.end / HASH_VAL)
                for i in range(start, end+1):
                    features_dict[chrom, i, element.strand].append(element)
            progress.update(1)
        return features_dict

    def _update_introns(self):
        for hash_val, features in self.features_dict.iteritems():
            for feature in features:
                if feature.featuretype == 'transcript':
                    for transcript_id in feature.attributes['transcript_id']:
                        exons = self.exons_dict[transcript_id]
                        transcript = self.transcripts_dict[transcript_id]
                        introns = self.find_introns(transcript, exons)
                        for intron in introns:
                            intron_feature = copy.deepcopy(feature)
                            intron_feature.start = intron['start']
                            intron_feature.end = intron['end']
                            intron_feature.featuretype = 'intron'
                            self.features_dict[hash_val].append(
                                intron_feature
                            )

    def _get_all_cds_dict(self):
        """
        For every cds-annotated transcript id (ENST), return a 
        dictionary containing the lowest and highest
        cds start and end vals for that transcript.

        :return cds_dict : defaultdict{transcript:{'start':START, 'end':END}} 
        """
        cds_dict = defaultdict(lambda: {'low': MAXVAL, 'hi': MINVAL})
        for cds_feature in self._db.features_of_type('CDS'):
            for transcript_id in cds_feature.attributes['transcript_id']:

                if cds_feature.start <= cds_dict[transcript_id]['low']:
                    cds_dict[transcript_id]['low'] = cds_feature.start
                if cds_feature.end >= cds_dict[transcript_id]['hi']:
                    cds_dict[transcript_id]['hi'] = cds_feature.end
        return cds_dict

    def _get_all_exons_dict(self):
        """
        :return:
        """
        exons_dict = defaultdict(list)
        for exon_feature in self._db.features_of_type('exon'):
            for transcript_id in exon_feature.attributes['transcript_id']:
                exons_dict[transcript_id].append(
                    {
                        'start': exon_feature.start,
                        'end': exon_feature.end
                    }
                )
        return exons_dict

    def _get_all_transcripts_dict(self):
        """

        :return:
        """
        transcripts_dict = defaultdict(list)
        for transcript_feature in self._db.features_of_type('transcript'):
            for transcript_id in transcript_feature.attributes['transcript_id']:
                transcripts_dict[transcript_id] = {
                    'start': transcript_feature.start,
                    'end': transcript_feature.end
                }
        return transcripts_dict

    def _infer_introns(self, transcript_id):
        """
        Returns a list of introns

        :param transcript_id:
        :return:
        """

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
                if self.cds_dict[transcript_id]['low'] > utr_feature.end:
                    five_prime_utr = True
                if self.cds_dict[transcript_id]['hi'] < utr_feature.start:
                    three_prime_utr = True
            elif utr_feature.strand == '-':
                if self.cds_dict[transcript_id]['low'] > utr_feature.end:
                    three_prime_utr = True
                if self.cds_dict[transcript_id]['hi'] < utr_feature.start:
                    five_prime_utr = True

        if five_prime_utr and three_prime_utr:
            return 'THREE_AND_FIVE_PRIME_UTR'
        elif five_prime_utr:
            return '5UTR'
            # return 'five_prime_utr'
        elif three_prime_utr:
            return '3UTR'
            # return 'three_prime_utr'
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
            gene_id = gene.attributes['gene_id'][0] if type(
                gene.attributes['gene_id']
            ) == list else gene.attributes['gene_id']
            try:
                gene_name_dict[gene_id] = gene.attributes['gene_name'][0]
            except KeyError:
                print(gene.attributes.keys())
                print("Warning. Key not found for {}".format(gene))
                return 1
        return gene_name_dict

    def find_introns(self, transcript, exons):
        positions = []
        introns = []
        for exon in exons:
            positions.append(exon['start'] - 1)
            positions.append(exon['end'] + 1)
        positions = sorted(positions)

        if positions[0] < transcript[
            'start']:  # there is no intron at the start of the feature
            positions.pop(0)
        else:
            positions.insert(transcript['start'])
        if positions[-1] > transcript['end']:
            positions.pop(-1)
        else:
            positions.append(transcript['end'])
        for i in range(0, len(positions) - 1, 2):
            introns.append({'start': positions[i], 'end': positions[i + 1]})
        return introns

    def get_all_overlapping_features_from_query(self, chrom, qstart, qend,
                                                strand):
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

    def annotate(self, interval):
        overlapping_features = self.get_all_overlapping_features_from_query(
            interval.chrom,
            interval.start,
            interval.end,
            interval.strand
        )
        to_append = ''  # full list of genes overlapping features
        transcript = defaultdict(list)
        for feature in overlapping_features:  # for each overlapping feature
            for transcript_id in feature.attributes[
                'transcript_id'
            ]:  # multiple genes can be associated with one feature
                transcript[transcript_id].append(
                    feature)  # append features to their respective genes
        for transcript, features in transcript.iteritems():
            for feature in features:
                # if 'protein_coding' not in feature.attributes['transcript_type']:
                #     if feature.featuretype == 'exon' or feature.featuretype == 'UTR':
                #         feature.featuretype = 'noncoding_exon'
                if feature.featuretype == 'UTR':
                    feature.featuretype = self._classify_utr(feature)
                to_append += "{}:{}:{}:{}:{}:".format(
                    transcript,
                    feature.start,
                    feature.end,
                    feature.strand,
                    feature.featuretype,
                )
                for t in feature.attributes['gene_id']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
                for t in feature.attributes['gene_name']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
                for t in feature.attributes['transcript_type']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + '|'
        return to_append[:-1]

def annotate(db_file, bed_file, out_file, chroms):
    """
    Given a bed6 file, return the file with an extra column containing
    '|' delimited gene annotations

    :param db_file:
    :param bed_file:
    :param out_file:
    :param chroms:
    :return:
    """
    annotator = Annotator(db_file, chroms)
    bed_tool = pybedtools.BedTool(bed_file)
    with open(out_file, 'w') as o:
        for interval in bed_tool:  # for each line in bed file
            annotation = annotator.annotate(interval)
            o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                interval.chrom, interval.start,
                interval.end, interval.name,
                interval.score, interval.strand,
                annotation  # we dont need last pipe
            ))


def annotate_with_genes(db_file, bed_file, out_file, chroms):
    """
    Given a bed6 file, return the file with an extra column containing
    '|' delimited gene annotations

    :param db_file:
    :param bed_file:
    :param out_file:
    :return:
    """
    annotator = Annotator(db_file, chroms)
    bed_tool = pybedtools.BedTool(bed_file)
    with open(out_file, 'w') as o:
        for interval in bed_tool: # for each line in bed file
            overlapping_features = annotator.get_all_overlapping_features_from_query(
                interval.chrom,
                interval.start,
                interval.end,
                interval.strand
            )
            to_append = '' # full list of genes overlapping features
            genes = defaultdict(list)
            for feature in overlapping_features: # for each overlapping feature
                for gene_id in feature.attributes['gene_id']: # multiple genes can be associated with one feature
                    genes[gene_id].append(feature)  # append features to their respective genes
            for gene in genes.keys():
                to_append = to_append + gene + '|'
            o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                interval.chrom, interval.start,
                interval.end, interval.name,
                interval.score, interval.strand,
                to_append[:-1] # we dont need last pipe
            ))