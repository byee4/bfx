__author__ = 'brian'

import pandas

def counts_to_rpkm(feature_counts_table):
    counts = feature_counts_table.ix[:,5:]
    lengths = feature_counts_table['Length']
    mapped_reads = counts.sum()
    return (counts * pow(10,9)).div(mapped_reads, axis=1).div(lengths, axis=0)

def counts_to_tpm(featureCountsTable):
    rpkm = counts_to_rpkm(featureCountsTable)
    tpm = rpkm.div(rpkm.sum())*pow(10,6)
    return tpm