__author__ = 'brian'

import pandas

def counts_to_rpkm(counts_table, skip_col=5):
    """
    simple function that converts a featureCounts pandas Dataframe 
    into an rpkm dataframe.
    
    :param counts_table: pandas.DataFrame()
        either a featureCounts table (first five cols contain non-count info,
        the rest contain raw counts) or a generic counts table (use skip_col=0
        in this case)
    :return rpkm: pandas.DataFrame 
    """
    counts = counts_table.ix[:, skip_col:]
    lengths = counts_table['Length']
    mapped_reads = counts.sum()
    return (counts * pow(10,9)).div(mapped_reads, axis=1).div(lengths, axis=0)


def counts_to_tpm(counts_table, skip_col=5):
    """
    simple function that converts a featureCounts pandas Dataframe
    into a TPM dataframe.
    
    :param counts_table: pandas.DataFrame() 
        either a featureCounts table (first five cols contain non-count info,
        the rest contain raw counts) or a generic counts table (use skip_col=0
        in this case)
    :return tpm: pandas.DataFrame
    """
    rpkm = counts_to_rpkm(counts_table, skip_col)
    tpm = rpkm.div(rpkm.sum())*pow(10,6)
    return tpm