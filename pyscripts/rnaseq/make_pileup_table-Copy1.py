'''
Created on Feb 2, 2017

@author: brianyee
'''
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import pysam
import seaborn as sns
import pybedtools
from collections import defaultdict
from Bio import SeqIO

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
    return df['A']+df['G']+df['C']+df['T']


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
    alphabet = ['A','C','G','T']
    for col, val in row.iteritems():
        if val > mx and col != row['ref'] and col in alphabet:
            mx = val
            mxcol = col
    return mxcol


def get_alt_cov(row):
    """get the total coverage of non-ref letters"""
    x = 0
    alphabet = ['A','C','T','G']
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


def ini(ref=''):
    alphabet = defaultdict(str)
    alphabet['A'] = 0
    alphabet['T'] = 0
    alphabet['C'] = 0
    alphabet['G'] = 0
    # alphabet['N'] = 0
    # alphabet['del'] = 0
    # alphabet['skip'] = 0
    alphabet['ref'] = ref
    return alphabet


def get_position_matrix(bam, chrom, start, maxlen, reffile, stepper='all'):
    total_reads = 0
    print('max len: {}'.format(maxlen))
    reference = pybedtools.BedTool.seq([chrom,0,maxlen],reffile)
    infile = pysam.AlignmentFile(bam, "rb", reference_filename=reffile)
    print('chrom: {}, start: {}, end: {}, len: {}'.format(chrom, 0, maxlen, len(reference)))
    count = start # running counter for each position added
    
    alphabet = ini()
    print("STARTING POS MATRIX")
    positions = [] # total length = maxlen - start
    offset = 0 # 
    max_offset = 0
    MAX_DEPTH = 10000000
    check = start
    for pileupcolumn in infile.pileup(chrom,start,maxlen,stepper=stepper,max_depth=MAX_DEPTH):
        if pileupcolumn.pos >= start: # if the read start occurs before the start of the chromosome or coordinate, skip:
            st = "" # foreach 
            if count >= (maxlen) or pileupcolumn.pos >= maxlen: # I think this works because of sorted reads?
                break
            """
            if there is no read coverage at the beginning positions
            """
            while count < pileupcolumn.pos:
                print(pileupcolumn.reference_pos, maxlen, len(reference)) # , reference[pileupcolumn.reference_pos].upper())
                alphabet = ini(reference[pileupcolumn.reference_pos].upper())
                
                positions.append(alphabet)
                alphabet = ini()
                count = count + 1

                # print('{}\t0'.format(count))
            # print str(pp.pos)+'\t'+str(pp.n)
            # print(len(pileupcolumn.pileups))
            for pileupread in pileupcolumn.pileups: # for each pileup read
                # offset = start - min(pileupread.alignment.get_reference_positions())
                # if offset > max_offset:
                #     max_offset = offset
                # if(pileupread.query_position >=offset):
                        # print(pileupcolumn.pos)
                total_reads = total_reads + 1
                if not pileupread.is_del and not pileupread.is_refskip:
                    alphabet[pileupread.alignment.query_sequence[pileupread.query_position].upper()] +=1
                    # st = st + pileupread.alignment.query_sequence[pileupread.query_position]
                elif pileupread.is_del:
                    alphabet['del'] +=1
                    # st = st + 'd'
                elif pileupread.is_refskip:
                    alphabet['skip'] +=1
                    # st = st + 's'
                else:
                    print("should never happen.")
                    # print(pileupread.is_head, pileupread.is_tail, pileupread.query_position, pileupread.indel, pileupread.alignment.query_sequence[pileupread.query_position])
                    # st = st + pileupread.alignment.query_sequence[pileupread.query_position]
            
                # print(st)
                # print("ADDING: {} at step: {}, at pos: {}".format(st, count, pileupcolumn.reference_pos))
            #alphabet['A'] = st.count('A')
            #alphabet['T'] = st.count('T')
            #alphabet['C'] = st.count('C')
            #alphabet['G'] = st.count('G')
            #alphabet['del'] = st.count('d')
            alphabet['ref'] = reference[pileupcolumn.reference_pos].upper()
            count = count + 1
            positions.append(alphabet)
            alphabet = ini()
    """
    If there are positions in the end without read coverage
    """
    while(count < maxlen):
        alphabet = ini(reference[count].upper())
        count = count + 1
        positions.append(alphabet)
    positions = pd.DataFrame(positions)
    del positions['N'] # we already call deletions
    return positions


def plot_mpileup(bam, chrom, start, end, ref, output_dir, exon_coords = [], stepper = 'all'):
    cols = sns.color_palette("hls", 6)
    filename = os.path.join(output_dir, os.path.basename(bam))
    if(1): #not os.path.exists(filename.replace('.bam','.tsv')): # if exists, don't make
        # print(os.path.basename(bam))
        separators = []
        df = get_position_matrix(bam, chrom, start, end, ref, stepper = stepper)
        df['alt_cov'] = df.apply(get_alt_cov,axis=1)
        df['ref_cov'] = df.apply(get_ref_cov,axis=1)
        
        """
        logic for concatenating exons
        """
        if len(exon_coords) > 0:
            dfx = pd.DataFrame()
            for exon in exon_coords:
                dfx = pd.concat([dfx,df[exon[0]:exon[1]]])
                separators.append(dfx.shape[0])
        else:
            dfx = df
        
        """
        draw the map
        """
        dfx.reset_index(inplace=True,drop=True)
        dfx.to_csv(filename.replace('.bam','.tsv'),sep='\t')
        
        fig, ax1 = plt.subplots(figsize=(100, 5))
        plt.bar(range(0,dfx.shape[0]), dfx['alt_cov'], color='white', width=1.0, label = 'alt')
        plt.bar(range(0,dfx.shape[0]), dfx['ref_cov'], bottom = dfx['alt_cov'], width=1.0, color='blue', label = 'ref')
        plt.bar(range(0,dfx.shape[0]), -dfx['del'], color=cols[3], width=1.0, label = 'del')
        for line in separators:
            plt.axvline(x=line)
        plt.legend()
        plt.savefig(filename.replace('.bam','.svg'), format='svg')
        plt.savefig(filename.replace('.bam','.png'), format='png')
        
        plt.close(fig)
        plt.gcf().clear()


def get_transcripts_from_genelist(genelist,gene_list = []):
    """
    given a genelist file, return the corresponding transcript
    """
    df = pd.read_table(genelist,names=['chrom','ensg','gene','transcript'])
    if(len(gene_list)==0):
        return list(df['transcript'])
    else:
        return list(df[df['gene'].isin(gene_list)])

######

def get_position_matrix_efficiently(bam, chrom, start, maxlen, reffile,
                                    stepper='all'):
    total_reads = 0

    infile = pysam.AlignmentFile(bam, "rb", reference_filename=reffile)

    count = start  # running counter for each position added

    alphabet = ini()
    # print("STARTING POS MATRIX")
    positions = []  # total length = maxlen - start
    offset = 0  #
    max_offset = 0
    MAX_DEPTH = 1000
    for pileupcolumn in infile.pileup(chrom, start, maxlen, stepper=stepper,
                                      max_depth=MAX_DEPTH):
        refbase = pybedtools.BedTool.seq(
            [chrom, pileupcolumn.reference_pos,
             pileupcolumn.reference_pos + 1],
            reffile
        )

        if pileupcolumn.pos >= start:  # if the read start occurs before the start of the chromosome or coordinate, skip:
            st = ""  # foreach
            if count >= (
            maxlen) or pileupcolumn.pos >= maxlen:  # I think this works because of sorted reads?
                break
            while count < pileupcolumn.pos:
                refbase = pybedtools.BedTool.seq(
                    [chrom, pileupcolumn.reference_pos,
                     pileupcolumn.reference_pos + 1],
                    reffile
                )
                alphabet = ini('NOCOV')

                positions.append(alphabet)
                alphabet = ini()
                count = count + 1
            for pileupread in pileupcolumn.pileups:  # for each pileup read
                if not pileupread.is_del and not pileupread.is_refskip:
                    alphabet[pileupread.alignment.query_sequence[
                        pileupread.query_position].upper()] += 1
            alphabet['ref'] = refbase.upper()
            count = count + 1
            positions.append(alphabet)
            alphabet = ini()
    while (count < maxlen):
        alphabet = ini('NOCOV')
        count = count + 1
        positions.append(alphabet)
    positions = pd.DataFrame(positions)
    return positions

def get_ref_num(row):
    ref_base = row['ref']
    return float(row[[ref_base]].to_string(index=False))

def get_alt_num(row):
    ref_num = get_ref_num(row)
    tot_num = row[['A','T','C','G']].sum()
    return tot_num - ref_num

def get_percent(row):
    if row['alt_num'] + row['ref_num'] == 0:
        return np.nan
    return row['alt_num'] / (row['alt_num'] + row['ref_num'])


def get_mut_percent(df):
    df['ref_num'] = df.apply(get_ref_num, axis=1)
    df['alt_num'] = df.apply(get_alt_num, axis=1)
    df['percent_mut'] = df.apply(get_percent, axis=1)
    return df['percent_mut']

def get_merged_table(bam, exons_file, ref, out_file):
    exons = pybedtools.BedTool(exons_file)
    e = exons[0]
    OFFSET = 70
    merged = get_position_matrix_efficiently(bam, str(e.chrom), e.start-OFFSET, e.start+OFFSET, ref, 'all')
    merged = get_mut_percent(merged)

    count = 0
    for e in exons:
        if e.end - e.start > OFFSET:
            df = get_position_matrix_efficiently(bam, str(e.chrom), e.start-OFFSET, e.start+OFFSET, ref, 'all')
            print(e.chrom, e.start-OFFSET, e.start+OFFSET)
            df = get_mut_percent(df)
            merged = pd.concat([merged, df], axis=1)
            count += 1
    merged.to_csv(out_file, sep='\t')

def main(argv=None): # IGNORE:C0111
    
    # Setup argument parser
    # USAGE: 
    # python plot_features_from_xintao_using_erics_manifest.py -o /projects/ps-yeolab3/bay001/maps/se/xintao/8-15-2016 -f -m /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS_submittedonly.txt -e se -r /projects/ps-yeolab3/bay001/maps/alt_splicing/8-5-2016/xintao-as-miso -s
    # manifest file is taken from here: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    
    parser.add_argument("-bam", "--bam", dest="bam",required=True)
    # parser.add_argument("-g", "--genelist", dest="genelist",required=True)
    # parser.add_argument("-o", "--outdir", dest="outdir",required=True)
    parser.add_argument("-o", "--outfile", dest="outfile", required=True)
    parser.add_argument("-s", "--stepper", dest="stepper", default="nofilter")
    parser.add_argument("-r", "--reffile", dest="reffile", required=True)
    # parser.add_argument("-c", "--chrom", dest="chrom",required=True)
    parser.add_argument("-e", "--exons", dest="exons", required=False)
    args = parser.parse_args()
    
    # output
    # outdir = args.outdir
    outfile = args.outfile
    bam = args.bam
    stepper = args.stepper
    reffile = args.reffile
    # chrom = args.chrom
    # genelist = args.genelist
    exons_file = args.exons
    """
    logger = logging.getLogger('map')
    logger.setLevel(logging.INFO)
    # ih = logging.FileHandler(os.path.join(outdir,'.map.txt'))
    # dh = logging.FileHandler(os.path.join(outdir,'.map.debug'))
    # eh = logging.FileHandler(os.path.join(outdir,'.map.err'))
    ih.setLevel(logging.INFO)
    dh.setLevel(logging.INFO)
    eh.setLevel(logging.ERROR)
    logger.addHandler(ih)
    logger.addHandler(dh)
    logger.addHandler(eh)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ih.setFormatter(formatter)
    dh.setFormatter(formatter)
    eh.setFormatter(formatter)
    """
    # chroms = get_transcripts_from_genelist(genelist)
    # logger.debug("Building the fasta dictionary")
    fasta_dict = build_fasta_dict(reffile)
    
    
    # start = 0
    # end = get_seq_len(fasta_dict,chrom)
    # logger.debug("Making pileup for {} - {}".format(os.path.basename(bam),chrom))
    get_merged_table(bam, exons_file, reffile, outfile)
    # plot_mpileup(bam, chrom, start, end, reffile, outdir, exon_coords = [], stepper = stepper)
if __name__ == "__main__":
    
    main()
