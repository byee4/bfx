'''
Created on July 7, 2017

Customized script for plotting the pileup of reads that map to repetitive elements.
At some point this will deprecate cut_sam, etc.

@author: brianyee
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pandas as pd
from collections import defaultdict
import os
import glob
import pysam
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools

### FROM CUT SAM ###
def get_element_dict(sam, elements):
    """
    From an element *RNU1, RNU2, etc., and an "eric" sam file, return a dictionary
    containing the reads corresponding to {element:[reads]}

    Let's ignore the _spliced transcripts for now...

    I expect the 'chrom' column to look like this:

    RNA28S||ENST00000607521.1
    chr10
    antisense_Alu||antisense_AluSx
    ...

    """
    chromosome_dict = defaultdict(list)
    with open(sam,'r') as f:
        for line in f:
            line = line.split('\t')
            for element in elements:
                if line[2].startswith("{}||".format(element)):
                    key = line[2]
                    line[2] = line[2].split('||')[1]
                    truncated_line = line[:12] # only grab the first 11 tabs
                    chromosome_dict[key].append('\t'.join(truncated_line))
    return chromosome_dict


def get_transcripts_from_genelist(genelist, gene_list = []):
    """
    given a genelist file, return the corresponding transcript.
    If gene_list is empty, return all transcripts.
    If gene_list contains a list of gene names, return transcripts
    only corresponding to the genes in this list.\

    see: pyscripts/data/genelist.txt
    """
    df = pd.read_table(genelist,names=['chrom','ensg','gene','transcript'])
    if(len(gene_list)==0):
        return list(df['transcript'])
    else:
        return list(df[df['gene'].isin(gene_list)])


def get_elements_from_genelist(genelist, gene_list = []):
    """
    'element' corresponds to the 'chrom' column of eric's rep element sam file.

    given a genelist file, return the corresponding transcript
    If gene_list is empty, return all transcripts.
    If gene_list contains a list of gene names, return transcripts
    only corresponding to the genes in this list.\

    see: pyscripts/data/genelist.txt
    """
    df = pd.read_table(genelist,names=['chrom','ensg','gene','transcript'])
    if(len(gene_list)==0):
        return list(df['chrom'])
    else:
        return list(df[df['gene'].isin(gene_list)])


def cut_sam(input_eric_file, output_file_prefix, output_sam_dir, gene_list):
    """
    Takes an eric-formatted sam-style file and splits each read into separate
    files specified with gene_list.

    :param input_eric_file: basestring
    :param output_sam_dir: basestring
    :param gene_list: basestring
    :return:
    """
    outfiles_created = []

    kept_transcripts = get_transcripts_from_genelist(gene_list) # all the x we want to keep
    elements = get_elements_from_genelist(gene_list) # all the e we want to keep

    chrdict = get_element_dict(input_eric_file, elements)
    for key, value in chrdict.iteritems():
        current_transcript = key.split('||')[1]
        out_file = os.path.join(
            output_sam_dir,
            os.path.basename(
                output_file_prefix + '.{}.cut.sam'.format(
                    current_transcript.replace('.','_')
                )
            )
        )
        if current_transcript in kept_transcripts:
            print('creating file: {}'.format(out_file))
            pd.DataFrame([line.split('\t') for line in value]).to_csv(
                out_file, sep='\t', index=None, header=None
            )
            outfiles_created.append(out_file)
    return outfiles_created

### OFFSET FUNCTIONS

def get_seq(allfastas, seq_header):
    """
    Given a sequence identifier, return the sequence.

    :param allfastas: basestring
    :param seq_header: basestring
    :return seq: basestring
    """
    handle = open(allfastas)

    for record in SeqIO.parse(handle, "fasta"):
        if seq_header == record.id:
            return str(record.seq)
    handle.close()

def create_offset_file(
        in_file,
        minor_name, minor_reference,
        major_name, major_reference,
        out_file
):
    """
    This function is basically a lazy way (.find()) that pushes (offsets)
    each read by re-mapping each individual element to a 'master' reference.
    This allows multiple rep elements to map to the same reference.

    :param in_file:
        input CUT sam file
    :param minor_name: basestring
        The name of the transcript contained within the major_reference
    :param minor_reference: basestring
        FASTA file of the transcript that contains all minor sequences
    :param major_name: basestring
        The name of the transcript that contains all minor sequences
    :param major_reference: basestring
        FASTA file of the transcript that contains the major sequence
    :param out_file:
        output CUT+OFFSET sam file
    :return:
    """

    offset = 0
    reference = get_seq(major_reference, major_name)

    with open(in_file, 'r') as f:

        offset = reference.find(
            get_seq(minor_reference, minor_name)
        )  # TODO: implement a proper aligner
        if offset != -1:
            o = open(out_file, 'w')
            for line in f:
                line = line.split('\t')
                line[2] = major_name # rename the CHROM column
                line[3] = str(int(line[3]) + offset) # modify the POS column
                o.write('\t'.join(line))
            o.close()
        else:
            print("WARNING: {} not found in {}. Skipping!".format(minor_name, major_name))

### FROM ADD SAM HEADER CONVERT BAM ###

def add_header_convert_to_bam(cut_sam_list, reference_fa):
    """
    Adds SAM headers and converts to bam

    :param cut_sam_list:
    :param reference_fa:
    :return:

    """
    bam_files = []
    for sam_file in cut_sam_list:
        bam_file = sam_file.replace('.sam','.bam')
        pysam.view(
            "-T", reference_fa, "-b", "-o",
            bam_file,
            sam_file, catch_stdout=False
        )
        bam_files.append(bam_file)

    return bam_files

### SORT BAM FILES ###

def sort_bam(bam_files):
    """
    Sorts all bam files in a directory given a sam prefix to search

    :param bam_files:
    :return:
    """
    sorted_files = []

    for bam_file in bam_files:
        sorted_file = bam_file.replace('.bam', '.sorted.bam')
        pysam.sort(
            "-o", sorted_file,
            bam_file
        )
        pysam.index(sorted_file)
        sorted_files.append(sorted_file)

    return sorted_files

### PLOT MPILEUP ###

def build_fasta_dict(fasta):
    """returns dictionary of fasta id:sequence"""
    fasta_dict = {}
    FastaFile = open(fasta, 'rU')
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.id
        seq = str(rec.seq)
        fasta_dict[name] = seq
    FastaFile.close()
    return fasta_dict


def get_seq_len(fasta_dict, key):
    return len(fasta_dict[key])


def get_cov(df):
    """returns total coverage as a series"""
    return df['A'] + df['G'] + df['C'] + df['T']


def get_ref(row):
    """returns the ref nucleotide at that row"""
    return row[row['ref']]


def get_del(df):
    """returns the number of deletions as a series"""
    return df['del']


def get_max(row):
    """
    not including the reference position,
    get the column(letter) with the most coverage.
    """
    mx = -1
    alphabet = ['A', 'C', 'G', 'T']
    for col, val in row.iteritems():
        if val > mx and col != row['ref'] and col in alphabet:
            mx = val
            mxcol = col
    return mxcol


def get_alt_cov(row):
    """get the total coverage of non-ref letters"""
    x = 0
    alphabet = ['A', 'C', 'T', 'G']
    for col, val in row.iteritems():
        if col != row['ref'] and col in alphabet:
            x = x + val
    return x


def get_del_cov(row):
    """get number of deletions at that position"""
    return row['del']


def get_ref_cov(row):
    """get ref coverage at that position"""
    return row[row['ref']]


def get_position_matrix(bam, chrom, start, maxlen, reffile, stepper='all'):
    total_reads = 0
    reference = pybedtools.BedTool.seq([chrom, 0, maxlen], reffile)
    infile = pysam.AlignmentFile(bam, "rb", reference_filename=reffile)
    count = start  # running counter for each position added
    alphabet = {}
    positions = []
    offset = 0
    max_offset = 0
    MAX_DEPTH = 10000000
    check = start
    for pileupcolumn in infile.pileup(chrom, start, maxlen, stepper=stepper,
                                      max_depth=MAX_DEPTH):
        if pileupcolumn.pos >= start:
            st = ""
            # print("pileuppos: {}".format(pileupcolumn.pos))
            # print("count: {}".format(count))
            if count >= (
            maxlen) or pileupcolumn.pos >= maxlen:  # I think this works because of sorted reads?
                break
            """
            if there is no read coverage at the beginning positions
            """
            while (count < pileupcolumn.pos):
                alphabet['A'] = 0
                alphabet['T'] = 0
                alphabet['C'] = 0
                alphabet['G'] = 0
                alphabet['del'] = 0
                alphabet['ref'] = reference[pileupcolumn.reference_pos].upper()
                # print(alphabet)
                positions.append(alphabet)
                alphabet = {}
                # print("ADDING COUNT")
                count = count + 1

                # print('{}\t0'.format(count))
            # print str(pp.pos)+'\t'+str(pp.n)
            # print(len(pileupcolumn.pileups))
            for pileupread in pileupcolumn.pileups:  # for each pileup read
                # offset = start - min(pileupread.alignment.get_reference_positions())
                # if offset > max_offset:
                #     max_offset = offset
                # if(pileupread.query_position >=offset):
                # print(pileupcolumn.pos)
                total_reads = total_reads + 1
                if not pileupread.is_del and not pileupread.indel and not pileupread.is_refskip:
                    st = st + pileupread.alignment.query_sequence[
                        pileupread.query_position]
                elif pileupread.is_del:
                    st = st + 'd'
                elif pileupread.is_refskip:
                    st = st + 's'
                else:
                    st = st + '-'

                    # print(st)
                    # print("ADDING: {} at step: {}, at pos: {}".format(st, count, pileupcolumn.reference_pos))
            alphabet['A'] = st.count('A')
            alphabet['T'] = st.count('T')
            alphabet['C'] = st.count('C')
            alphabet['G'] = st.count('G')
            alphabet['del'] = st.count('d')
            alphabet['ref'] = reference[pileupcolumn.reference_pos].upper()
            count = count + 1
            # print(alphabet)

            positions.append(alphabet)
            alphabet = {}
            # print('{} '.format(count)),
    """
    If there are positions in the end without read coverage
    """
    while (count < maxlen):
        # count = count + 1
        alphabet['A'] = 0
        alphabet['T'] = 0
        alphabet['C'] = 0
        alphabet['G'] = 0
        alphabet['del'] = 0
        alphabet['ref'] = reference[count].upper()
        # print(alphabet)
        count = count + 1
        positions.append(alphabet)
    # print(start, maxlen, len(positions), max_offset, check-start)
    # print(total_reads)
    return pd.DataFrame(positions)


def plot_mpileup(bam, chrom, start, end, ref, output_dir, exon_coords=[],
                 stepper='all'):
    cols = sns.color_palette("hls", 6)
    filename = os.path.join(output_dir, os.path.basename(bam))
    if not os.path.exists(
            filename.replace('.bam', '.tsv')):  # if exists, don't make
        # print(os.path.basename(bam))
        separators = []
        df = get_position_matrix(bam, chrom, start, end, ref, stepper=stepper)
        df['alt_cov'] = df.apply(get_alt_cov, axis=1)
        df['ref_cov'] = df.apply(get_ref_cov, axis=1)

        """
        logic for concatenating exons
        """
        if len(exon_coords) > 0:
            dfx = pd.DataFrame()
            for exon in exon_coords:
                dfx = pd.concat([dfx, df[exon[0]:exon[1]]])
                separators.append(dfx.shape[0])
        else:
            dfx = df

        """
        draw the map
        """
        dfx.reset_index(inplace=True, drop=True)
        fig, ax1 = plt.subplots(figsize=(100, 5))
        plt.bar(range(0, dfx.shape[0]), dfx['alt_cov'], color='white',
                width=1.0, label='alt')
        plt.bar(range(0, dfx.shape[0]), dfx['ref_cov'], bottom=dfx['alt_cov'],
                width=1.0, color='blue', label='ref')
        plt.bar(range(0, dfx.shape[0]), -dfx['del'], color=cols[3], width=1.0,
                label='del')
        for line in separators:
            plt.axvline(x=line)
        plt.legend()
        plt.savefig(filename.replace('.bam', '.svg'), format='svg')
        plt.savefig(filename.replace('.bam', '.png'), format='png')
        dfx.to_csv(filename.replace('.bam', '.tsv'), sep='\t')
        plt.close(fig)
        plt.gcf().clear()

def main(argv=None): # IGNORE:C0111
    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("-sam", "--sam", dest="sam", required=True)
    parser.add_argument("-g", "--genelist", dest="genelist", required=True)
    parser.add_argument("-o", "--outdir", dest="outdir", required=True)
    parser.add_argument("-r", "--reference", dest="reference", required=True)

    args = parser.parse_args()

    # output
    outdir = args.outdir
    sam = args.sam
    gene_list = args.genelist
    ref = args.reference

    prefix = os.path.splitext(os.path.basename(sam))[0]
    print('prefix: {}'.format(prefix))
    cut_sams = cut_sam(sam, prefix, outdir, gene_list)
    bam_files = add_header_convert_to_bam(cut_sams, ref)
    sorted_bams = sort_bam(bam_files)

if __name__ == "__main__":
    main()
