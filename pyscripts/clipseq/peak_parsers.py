import pybedtools
import os
import functools
import glob
import pandas as pd

READSBYLOC_COMBINED_CSV_L2FCWITHPVAL_ENR_HEADERS_1IP1IN = [ # deprecated
    'ENSG','CDS','CDS-pvalue','5utr','5utr-pvalue','3utr','3utr-pvalue',
    '5_and_3_utr','5_and_3_utr-pvalue','intron','intron-pvalue',
    'intergenic','intergenic-pvalue','noncoding_exon','noncoding_exon-pvalue',
    'noncoding_intron','noncoding_intron-pvalue'
]
READS_BY_LOC_HEADER = [ # deprecated
    'ENSG', 'CDS', '5utr', '3utr', '5utr|3utr', 'intron',
    'intergenic', 'noncoding_exon', 'noncoding_intron'
]

ANNOTATED_BED_HEADERS = [
    'chrom', 'start', 'end', 'pv', 'fc', 'strand', 'annotation', 'gene'
]

REGIONS = ['noncoding_exon', '3utr', '5utr', 'intron', 'noncoding_intron',
           'CDS', 'intergenic', '5utr_and_3utr']

import glob

def split_single_cols(df, col, sep='|'):
    """
    Splits a df['col'] into two separated by 'sep'
    """
    df["{} l2fc".format(col.split(sep)[1])], \
    df["{} l10p".format(col.split(sep)[1])] = zip(
        *df[col].map(lambda x: x.split(sep))
    )
    return df


def split_l2fcwithpval_enr(df, discard = True):
    """
    Splits a dataframe into its l2fc and log10 pvalue

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


def scatter_matrix(ip_l2fc, inp_reads_by_loc):
    """
    Inner joins the ip l2fc and input reads by loc files

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
    """ Doesn't really work if the header order changes
    plot_y = pd.read_table(
        ip_l2fc,
        sep='[\t\|]+',
        engine='python',
        index_col=0,
        skiprows=1,
        names=READSBYLOC_COMBINED_CSV_L2FCWITHPVAL_ENR_HEADERS_1IP1IN
    )"""

    plot_y = pd.read_table(
        ip_l2fc,
        index_col=0,
    )
    plot_y = split_l2fcwithpval_enr(plot_y)

    x = pd.merge(plot_x, plot_y, how='inner', left_index=True,
                 right_index=True)
    return x


def filter_input_norm_as_df(file_name, l2fc, pval, out_file=None):
    """
    Convenience method for me.

    Filters an "input norm"-formatted file given
    log2 fold change and log10 pvalue thresholds.
    See data/input_norm_bed.bed file for an example.

    Parameters
    ----------
    file_name : basename
    l2fc : float
    pval : float
    out_file : basename

    Returns
    -------
    filtered : pybedtools.BedTool()

    """
    return filter_input_norm(file_name, l2fc, pval, out_file).to_dataframe()


def filter_input_norm(file_name, l2fc, pval, out_file=None):
    """
    Filters an "input norm"-formatted file given
    log2 fold change and log10 pvalue thresholds.
    See data/input_norm_bed.bed file for an example.

    Parameters
    ----------
    file_name : basename
    l2fc : float
    pval : float
    out_file : basename

    Returns
    -------
    filtered : pybedtools.BedTool()

    """
    try:
        bedtool = pybedtools.BedTool(file_name)

        filter_data_inst = functools.partial(filter_data, l2fc=l2fc, pval=pval)

        bedtool = bedtool.filter(filter_data_inst)

        if out_file is not None:
            bedtool.saveas(out_file)

        return bedtool

    except Exception as e:
        return 1


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


def bed6_to_bed8(interval):
    """
    Basically appends the start/stop fields to 'thickStart/thickStop' fields
    Turns BED6 into BED8 (formerly called: make_clipper_ish)

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


def return_region(row):
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


def get_counts(wd, out_dir, l2fc, l10p, suffix='.annotated'):
    """
    Returns the number of peak counts for all regions
    annotated by eric's pipeline.

    Parameters
    ----------
    wd : string
        directory where the input_norm output is kept
        (where to look for .annotated files)
    """
    samples = {}

    for f in glob.glob(os.path.join(wd, '*{}'.format(suffix))):
        basename = os.path.basename(f)
        out_file = os.path.join(
            out_dir,
            basename.replace(
                '{}'.format(suffix),
                '{}.filtered-p{}-f{}'.format(
                    suffix, l10p, l2fc
                )
            ),
        )
        df = filter_input_norm_as_df(f, out_file, l2fc, l10p)
        df.columns = ANNOTATED_BED_HEADERS


        samples[basename] = {}
        df['region'] = df.apply(return_region, axis=1)
        for key, value in df['region'].value_counts().iteritems():
            samples[basename][key] = value
        for region in REGIONS:
            if region not in samples[basename]:
                samples[basename][region] = 0
    return pd.DataFrame(samples)