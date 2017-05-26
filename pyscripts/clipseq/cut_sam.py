'''
Created on Feb 2, 2017

Customized script for parsing eric's repetitive element mapped sam files
and turning them into proper sam files. Eric's files have weird chromosome fields
and other extra fields which embed repetitive family-mapped reads and other junk
but need to be SAM files in order to be properly used with other third party 
software.

@author: brianyee
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import os
import glob
import logging


def get_element_dict(sam, elements):
    """
    From an element *RNU1, RNU2, etc., and an "eric" sam file, return a dictionary
    Let's ignore the _spliced transcripts for now...
    """
    chromosome_dict = {}
    with open(sam,'r') as f:
        for line in f:
            line = line.split('\t')
            for element in elements:
                if(line[2].startswith("{}||".format(element))):
                    key = line[2]
                    if not key in chromosome_dict:
                        chromosome_dict[key] = []
                    line[2] = line[2].split('||')[1]
                    truncated_line = line[:12] # only grab the first 11 tabs
                    chromosome_dict[key].append('\t'.join(truncated_line))
    return chromosome_dict


def get_transcripts_from_genelist(genelist,gene_list = []):
    """
    given a genelist file, return the corresponding transcript
    """
    df = pd.read_table(genelist,names=['chrom','ensg','gene','transcript'])
    if(len(gene_list)==0):
        return list(df['transcript'])
    else:
        return list(df[df['gene'].isin(gene_list)])


def get_elements_from_genelist(genelist,gene_list = []):
    """
    given a genelist file, return the corresponding transcript
    """
    df = pd.read_table(genelist,names=['chrom','ensg','gene','transcript'])
    if(len(gene_list)==0):
        return list(df['chrom'])
    else:
        return list(df[df['gene'].isin(gene_list)])


def main(argv=None): # IGNORE:C0111
    
    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    
    parser.add_argument("-sam", "--sam", dest="sam",required=True)
    parser.add_argument("-g", "--genelist", dest="genelist",required=True)
    parser.add_argument("-o", "--outdir", dest="outdir",required=True)
    
    args = parser.parse_args()
    
    # output
    outdir = args.outdir
    sam = args.sam
    genelist = args.genelist
    
    logger = logging.getLogger('cut')
    logger.setLevel(logging.INFO)
    ih = logging.FileHandler(os.path.join(outdir,'cut.txt'))
    dh = logging.FileHandler(os.path.join(outdir,'cut.debug'))
    eh = logging.FileHandler(os.path.join(outdir,'cut.err'))
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
    
    kept_transcripts = get_transcripts_from_genelist(genelist)    
    elements = get_elements_from_genelist(genelist)
    
    try:
        o = glob.glob(outdir + os.path.basename(sam).replace('.sam','') + "*.cut.sam")
        if(len(o)>0): # if we've made at least one...
            for i in range(0,len(o)):
                logger.info('SUCCESS: {} exists'.format(o[i]))
                print(','),
        else:
            chrdict = get_element_dict(sam,elements)
            for key, value in chrdict.iteritems():
                current_transcript = key.split('||')[1]
                outfile = os.path.join(outdir,os.path.basename(sam.replace('.sam','.{}.cut.sam'.format(current_transcript.replace('.','_')))))
                if current_transcript in kept_transcripts:
                    pd.DataFrame([line.split('\t') for line in value]).to_csv(outfile, # replace . with _, not sure if this is the best idea but 
                    sep='\t',index=None,header=None)
                    logger.info('SUCCESS: {}'.format(outfile))
                    print('.'),
    except Exception as e:
        print(':'),
        logger.error('FAIL: {}'.format(sam))
if __name__ == "__main__":
    
    main()
