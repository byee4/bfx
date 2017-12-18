import pybedtools
import pandas as pd
import numpy as np
import cPickle as pickle
import os
from collections import defaultdict

def split_single_cols(df, col, sep='|'):
    """
    Splits a df['col'] into two separated by 'sep'
    ie. -0.9201|0.00000 -> -0.9201  0.00000
    """
    df["{} l2fc".format(col.split(sep)[1])], \
    df["{} l10p".format(col.split(sep)[1])] = zip(
        *df[col].map(lambda x: x.split(sep))
    )
    return df


def split_l2fcwithpval_enr(df, discard = True):
    """
    Splits a dataframe into its l2fc and log10 pvalue
    ie. -0.9201|0.00000 -> -0.9201  0.00000
    
    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of l2fcwithpval_enr file
    discard : bool
        if True, discard the original column.
        if False, keep the original column.

    Returns
    -------
    df : pandas.DataFrame
    """

    for col in df.columns:
        df = split_single_cols(df, col)
        if discard:
            del df[col]
    return df


def read_l2fcwithpval_enr(fn):
    """
    Reads in a *l2fcwithpval_enr.csv file and returns a dataframe.

    :param fn: 
    :return: 
    """
    df = pd.read_table(
        fn,
        index_col=0,
    )
    df = split_l2fcwithpval_enr(df)
    df = df.replace('NaN', np.nan)
    df = df.apply(pd.to_numeric)
    return df

def scatter_matrix(ip_l2fc, inp_reads_by_loc):
    """
    Inner joins the ip l2fc and input reads by loc files
    and builds a dataframe containing: input reads by location, 
    ip 
    Parameters
    ----------
    ip_l2fc : basestring
        "IP*_ReadsByLoc_combined.csv.l2fcwithpval_enr.csv"
    inp_reads_by_loc : basestring
        "INPUT*reads_by_loc.csv"

    Returns
    -------
    x : pandas.DataFrame
        intersection of fold-enrichment and input reads covering each gene.
    """

    plot_x = pd.read_table(
        inp_reads_by_loc,
        index_col=0,
    )

    plot_y = read_l2fcwithpval_enr(ip_l2fc)

    x = pd.merge(plot_x, plot_y, how='inner', left_index=True,
                 right_index=True)
    return x

def filter_input_norm(file_name, l10p, l2fc, as_bedtool=False):
    """
    Filters an "input norm"-formatted file given
    log2 fold change and log10 pvalue thresholds.
    See data/input_norm_bed.bed file for an example.
    
    :param file_name: basestring
        filename of the input normalized file.
    :param l2fc: float
        log2 fold change cutoff
    :param l10p: float
        -log10 pvalue cutoff
    :param col_names: list
        column names of the bed or bedlike file.
    :param as_bedtool: bool
        if True, return as bedtool. 
        if False, return dataframe.
    :return: 
    """
    try:
        col_names = ['chrom','start','end','l10p','l2fc','strand']
        df = pd.read_table(file_name, names=col_names)
        dff = filter_df(df, l10p, l2fc)
        print("before: {}, after: {}".format(df.shape[0], dff.shape[0]))
        if as_bedtool:
            return pybedtools.BedTool.from_dataframe(dff)
        return dff

    except Exception as e:
        print(e)
        return 1

def filter_df(df, l10p, l2fc):
    """
    Returns a dataframe filtered using l10p and l2fc values.
    Assumes the columns are aptly named 'l2fc' and 'l10p'.
    
    :param df: 
    :param l2fc: 
    :param l10p: 
    :return: 
    """
    return df[(df['l2fc']>=float(l2fc)) & (df['l10p']>=l10p)]

def return_region_eric(row):
    """
    Given a row of a inputnormed bedfile, return region
    Row must be in the same format as a line in Eric's
    *.annotated file.

    """
    try:
        if row['annotation'] == 'intergenic':
            return 'intergenic'
        region = row['annotation'].split('|')[0]

        return region
    except Exception as e:
        print(e)
        print(row)

def read_annotated_file(fn, headers=None, src='brian'):
    """
    Reads in an annotated bedfile from either eric or me.
    Ensures that 'region' column is defined and set.
    Returns dataframe
    
    :param fn: basestring
        filename
    :param headers: list
        if src isn't from eric or me, you have to 
        explicitly set the header columns.
    :param src: basestring
        either 'brian' or 'eric' depending on whose script you use.
    :return df: pandas.DataFrame
    """
    if src == 'brian':
        headers=[
            'chrom','start','end','l10p','l2fc','strand','geneid',
            'genename','region','alloverlapping'
        ]
        df = pd.read_table(fn, names=headers)
    elif src == 'eric':
        headers=[
            'chrom','start','end','l10p','l2fc',
            'strand','annotation','geneid'
        ]
        df = pd.read_table(fn, names=headers)
        df['region'] = df.apply(return_region_eric, axis=1)
    else:
        assert 'region' in headers
        df = pd.read_table(fn, names=headers)
    return df

def get_region_counts(fn, headers, src):
    """
    Returns a dataframe of REGIONS and the number of peaks
    associated with each region.
    
    :param fn: 
    :return: 
    """
    return pd.DataFrame(
        read_annotated_file(fn, headers, src)['region'].value_counts()
    )

def get_counts_df(dfs):
    """
    Returns the counts of each region in the annotation file.
    
    :param dfs: dict
        dict of dataframes and their associated region value counts.
        Keys = filenames, values = dataframes belonging to that filename.
    :param headers: list
        optional if either 'brian' or 'eric' are src. Otherwise, 
        this MUST contain a 'region' category. 
    :param src: basestring
        either 'brian' or 'eric' to denote the annotation structure.
    :return merged: pandas.DataFrame
        table containing the sum of counts for each region in the annotated file.
    """
    columns = dfs.keys()
    assert len(columns) > 0
    annotated_df = dfs[columns[0]]
    merged = annotated_df
    # merged.columns = [os.path.basename(annotated_file)]
    for key in columns[1:]:
        merged = pd.merge(merged, dfs[key], how='left', left_index=True,
                          right_index=True)
    return merged

def get_counts(fns, headers=None, src='brian', basename=True):
    dfs = {}
    for fn in fns:
        if basename:
            key = os.path.basename(fn)
        else:
            key = fn
        dfs[key] = get_region_counts(fn, headers, src)
        dfs[key].columns = [key]
    return get_counts_df(dfs)


def remove_peaks(wt_peaks, ko_peaks, as_bedtool=False):
    """
    Given a bed file, remove any peak contained with a 'peaks_to_remove'
    file.

    Parameters
    ----------
    wt_peaks : pybedtools.BedTool
    ko_peaks : pybedtools.BedTool
    as_bedtool : bool

    Returns
    -------

    """
    # ko = pybedtools.BedTool(ko_peaks)
    # wt = pybedtools.BedTool(wt_peaks)


    bedtool = wt_peaks.intersect(ko_peaks, v=True)
    print('number of peaks before: {}, after: {}'.format(len(wt_peaks), len(bedtool)))

    if as_bedtool:
        return bedtool
    df = bedtool.to_dataframe()
    return df


def read_parsed(fn):
    """
    Reads Eric's parsed file from the repetitive element pipeline.

    Parameters
    ----------
    fn : basestring
        the *.parsed file

    Returns
    -------
    total_df : pandas.DataFrame
        dataframe of total reads per unique/repetitive element family.
    element_df : pandas.DataFrame
        dataframe of unique repetitive/unique elements that each unique read
        mapped to.
    """
    df = pd.read_table(fn, comment='#', names=[
        'total_or_element','element','read_num',
        'percentage','annotation','gene'
    ])
    total_df = df[df['total_or_element']=='TOTAL'][
        ['element','read_num','percentage']
    ]
    element_df = df[df['total_or_element']=='ELEMENT'][
        ['element','read_num','percentage']
    ]
    return total_df, element_df


# LEGACY functions to handle some of the old gabe and eric stuff #

def read_kmer_enrichment_from_pickle(
        pickle_file, region='all', k=6, col_name='zscore delta'
):
    """
    Reads in a pickle file from gabe's clip_analysis script and returns a
    dataframe containing kmers and their enriched z-scores

    :param pickle_file: basestring
        pickle filename output from gabe's clip_analysis script.
    :param region: basestring
        one of:
        'all', 'three_prime_utrs', 'five_prime_utrs', 'distintron500',
        'cds', 'proxintron500'
    :return df: pandas.DataFrame
    """
    loaded = pickle.load(open(pickle_file, 'rb'))
    df = pd.DataFrame(loaded['kmer_results'][region][k]).T
    df.columns = ['fg', 'bg', col_name]
    return df[[col_name]]

def bed6_to_bed8(interval):
    """
    Basically appends the start/stop fields to 'thickStart/thickStop' fields
    Turns BED6 into BED8 (formerly called: make_clipper_ish)
    (Helps with plot_clip_analysis_figures.py in CLIPper). 
    
    Parameters
    ----------
    interval : pybedtools.Interval

    Returns
    -------

    """
    interval.name = interval[7]
    interval[6] = interval.start
    interval[7] = interval.stop

    return interval


def filter_data(interval, l2fc, pval):
    """
    col4 is -log10 p-val
    col5 is -log2 fold enrichment

    Expects the standard input norm file format.

    Parameters
    ----------
    interval : pybedtools.Interval
    l2fc : float
    pval : float

    Returns
    -------

    """

    return (float(interval[4]) >= pval) and (float(interval[3]) >= l2fc)


def split_annotated_file_into_region_beds_and_save(
    annotated, headers=None, src='brian',
    split_beds_directory=os.getcwd()
):
    if not os.path.exists(split_beds_directory):
        os.mkdir(split_beds_directory)

    df = read_annotated_file(annotated, headers, src)
    regions = set(df['region'])
    for region in regions:
        region_df = df[df['region']==region][[
            'chrom','start','end','l10p','l2fc','strand'
        ]]
        region_df.to_csv(
            os.path.join(
                split_beds_directory,
                os.path.basename(annotated) + ".{}.bed".format(
                    region
                )
            ),
            sep='\t',
            header=False,
            index=False
        )

### This is modified from Emily's code ###

def get_shuffled_region(bed_file_tool, bed_file_to_shuffle_within):
    """
    Makes a shuffled file
    :param bed_file: Bed file of peaks to shuffle
    :param bed_file_to_shuffle_within: Bed file of REGIONS to include in the shuffling.
        This is typicalls a file that matches the genic region of peaks (exons, 3'UTRs, etc.)
    :return: a shuffled bedtool object
    """

    # bed_file_tool = pybedtools.BedTool(bed_file)
    shuffled = bed_file_tool.shuffle(genome="hg19",
                                     incl=bed_file_to_shuffle_within,
                                     noOverlapping=True)
    return shuffled