'''
This is just the python runner for what was originally in a notebook.
Will need to clean this up later (if there's time.. haha)

@author: brianyee
'''

import glob
import manifest_helpers as m
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from collections import defaultdict
from tqdm import trange


__version__ = '0.0.1'


def get_dpsi_all(df):
    return get_dpsi_series(df, 1, 1, 1, -1)


def get_junction_region_from_row(row, event = 'se'):
    """
    Returns a stringified junction region for a row in rmats

    Parameters
    ----------
    row : pandas.core.series.Series

    Returns
    -------
    junction : basestring

    """
    if event == 'se':
        return '{}:{}-{}:{}'.format(
            row['chr'],
            row['upstreamEE'],
            row['downstreamES'],
            row['strand']
        )
    else:
        print('not implemented')
        return ''


def get_junction_region(df, event = 'se'):
    """
    Returns stringified junction regions for the whole rmats dataframe.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
    event : basestring
        currently only 'se' is supported.

    Returns
    -------
    df : pandas.core.frame.DataFrame

    """
    return df.apply(get_junction_region_from_row, axis=1, args=[event])


def get_avg_dpsi_for_all_junctions(df, default_region='junction'):
    """
    Returns the average dpsi from all junctions as a DataFrame

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
    default_region : basestring
        the region that describes the junctions in an rmats file
        (will do a groupby on this region and take the mean)

    Returns
    -------
    mean_dpsi : pandas.core.frame.DataFrame

    """
    return pd.DataFrame(df.groupby(default_region)['IncLevelDifference'].mean())


def filter_rmats_df(df, pv=1, fdr=1, inc_max=1, inc_min=-1, index=True):
    df = df[
        (df['PValue'] < pv) & (df['FDR'] < fdr) & \
        (df['IncLevelDifference'] <= inc_max) & \
        (df['IncLevelDifference'] >= inc_min)
        ]
    if index:
        df['index'] = df.apply(event_to_string, axis=1)
        df.set_index('index', inplace=True)
    return df

def get_dpsi_series(df, pv=1, fdr=1, inc_max=1, inc_min=-1, index=True):
    """

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        dataframe of the rmats JunctionCountOnly.txt file
    pv : float
        only values (<) pv will be returned
    fdr : float
        only values (<) fdr will be returned
    inc_max : float
        max (<=) incleveldifference to return
    inc_min : float
        min (>=) incleveldifference to return
    index : boolean
        determines whether or not we want the index preserved
    Returns
    -------
    dpsi_series : pandas.Series
        series list of all dpsi in the file
    """

    df = filter_rmats_df(df, pv, fdr, inc_max, inc_min, index)
    return df['IncLevelDifference']

def event_to_string(row, event = 'se'):
    """
    Returns a stringified version of an rmats row event.
    Events are chromosome sort ordered.

    Parameters
    ----------
    row
    event

    Returns
    -------

    """
    if event == 'se':
        return '{}@{}-{}@{}-{}@{}-{}@{}'.format(
            row['chr'],
            row['upstreamES'],
            row['upstreamEE'],
            row['exonStart_0base'],
            row['exonEnd'],
            row['downstreamES'],
            row['downstreamEE'],
            row['strand']
        )
    else:
        print("unimplemented")
        return ''

def get_events(df, event = 'se'):
    """
    Returns the event start and end for each described exon.

    Parameters
    ----------
    df

    Returns
    -------

    """
    if event == 'se':
        return df.apply(event_to_string, axis=1)

def read_rmats_df(f):
    """
    Returns the dataframe of f if exists, otherwise returns None.
    """
    if not os.path.exists(f):
        return None
    else:
        df = pd.read_table(f, index_col=0)
        return df


def get_avg_inclevels(row, level='IncLevel2'):
    """
    Returns the average of the two replicates for a given IncLevel column.
    Returns np.nan if both reps have NA values
    Returns the rep_x incLevel if rep_y incLevel is NA.

    Parameters
    ----------
    row :
    level

    Returns
    -------

    """
    tot = 0
    c = 0
    levels = row[level].split(',')

    for l in levels:
        if not l == 'NA':
            tot = tot + float(l)
            c += 1
    try:
        return tot / c
    except ZeroDivisionError as e:
        # number of valid samples is 0
        return np.nan


def get_jc_from_row(row, sample='SJC_SAMPLE_1'):
    """
    Returns the junction counts for the row for a given column (sample).

    Parameters
    ----------
    row : pandas.core.series.Series
    sample : basestring
        either one of: 'SJC_SAMPLE_1', 'SJC_SAMPLE_2', 'IJC_SAMPLE_1', 'IJC_SAMPLE_2'

    Returns
    -------
    r1, r2 : tuple
        rep1 and rep2 Skipped Junction Counts
    """
    r1, r2 = [int(i) for i in row[sample].split(',')]
    return r1, r2


def get_avg_jc_from_row(row, sample='SJC_SAMPLE_1'):
    """
    Returns the average junction counts for the row for a given column (sample)

    Parameters
    ----------
    row : pandas.core.series.Series
    sample : basestring
        either one of: 'SJC_SAMPLE_1', 'SJC_SAMPLE_2', 'IJC_SAMPLE_1', 'IJC_SAMPLE_2'

    Returns
    -------

    """
    r1, r2 = get_jc_from_row(row, sample)
    return (r1+r2)/2.0


def get_jc(df, sample='SJC_SAMPLE_1'):
    """
    Returns all jc tuples (r1, r2) for each row in the dataframe.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
    sample : basestring
        either one of: 'SJC_SAMPLE_1', 'SJC_SAMPLE_2', 'IJC_SAMPLE_1', 'IJC_SAMPLE_2'
    Returns
    -------

    """
    return df.apply(get_jc_from_row, axis=1, args=[sample])


def get_avg_jc(df, sample='SJC_SAMPLE_1'):
    """
    Returns the average jc for each row in the dataframe.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
    sample : basestring
        either one of: 'SJC_SAMPLE_1', 'SJC_SAMPLE_2', 'IJC_SAMPLE_1', 'IJC_SAMPLE_2'

    Returns
    -------
    junction_counts : pandas.core.series.Series

    """
    return df.apply(get_avg_jc_from_row, axis=1, args=[sample])


def get_inclevels(df, level='IncLevel1', rep=0):
    """
    Given a level and a rep number, return a series of inclevels from the df.
    NOTE: returns all NA's with np.nan !
    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        dataframe from rmats file.
    level : basestring
        either one of: 'IncLevel1' or 'IncLevel2'
    rep : int
        either one of: 0 or 1
    Returns
    -------
    inclevels : pandas.core.series.Series
    """
    return df.apply(get_inclevels_from_row, axis=1, args=([level, rep]))


def get_inclevels_from_row(row, level='IncLevel1', rep=0):
    """
    For a given row, an column (IncLevel1 or 2), and a replicate,
    return the incLevel.
    NOTE: returns all NA's with np.nan !
    """
    inc_level = row[level].split(',')[rep]
    if inc_level != 'NA':
        return float(inc_level)
    else:
        return np.nan


def split_and_append_all_csv(df):
    """
    Splits up all of the comma separated values in the table for easier
    parsing downstream.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        dataframe from rmats file.
    Returns
    -------
    df : pandas.core.frame.DataFrame
    """
    df['IJC_SAMPLE_1_avg'] = get_avg_jc(df, 'IJC_SAMPLE_1')
    df['IJC_SAMPLE_2_avg'] = get_avg_jc(df, 'IJC_SAMPLE_2')
    df['SJC_SAMPLE_1_avg'] = get_avg_jc(df, 'SJC_SAMPLE_1')
    df['SJC_SAMPLE_2_avg'] = get_avg_jc(df, 'SJC_SAMPLE_2')
    df['IncLevel1_rep1'] = get_inclevels(df, 'IncLevel1', 0)
    df['IncLevel1_rep2'] = get_inclevels(df, 'IncLevel1', 1)
    df['IncLevel2_rep1'] = get_inclevels(df, 'IncLevel2', 0)
    df['IncLevel2_rep2'] = get_inclevels(df, 'IncLevel2', 1)
    return df



def trim(dx, ext_label='all'):
    """
    calculates the inclevel for all reps and conditions in a single dataframe.
    Also labels the entire dataframe by the ext_label, so that multiple
    rmats files can be concatenated together without losing where each value came from.
    """
    inc1_rep1 = pd.DataFrame(
        dx.apply(get_inclevels_from_row, axis=1, args=('IncLevel1', 0))
    )
    inc1_rep1['label'] = 'IncLevel1 (KD??)'
    inc1_rep2 = pd.DataFrame(
        dx.apply(get_inclevels_from_row, axis=1, args=('IncLevel1', 1))
    )
    inc1_rep2['label'] = 'IncLevel1 (KD??)'
    inc2_rep1 = pd.DataFrame(
        dx.apply(get_inclevels_from_row, axis=1, args=('IncLevel2', 0))
    )
    inc2_rep1['label'] = 'IncLevel2 (Control??)'
    inc2_rep2 = pd.DataFrame(
        dx.apply(get_inclevels_from_row, axis=1, args=('IncLevel2', 1))
    )
    inc2_rep2['label'] = 'IncLevel2 (Control??)'
    inc = pd.concat([inc1_rep1, inc1_rep2, inc2_rep1, inc2_rep2])
    inc[0] = inc[0].replace(-1, np.nan)  # TODO: remove.
    inc['subset'] = ext_label
    return inc


def get_all_rmats(f, pos_suffix='.positive.nr.txt', neg_suffix='.negative.nr.txt'):
    """
    Given an rmats file (f) and using my special suffixes to distinguish
    between positive and negatively significant changing events, returns a list
    of 1) the original rmats file, 2) the significantly positively changing
    events, and 3) the significantly negatively changing events.
    """
    return [
        read_rmats_df(f),
        read_rmats_df(f.replace('.txt', pos_suffix)),
        read_rmats_df(f.replace('.txt', neg_suffix))
    ]


def add_labels_and_concat(df1, df2, df3):
    """
    Takes 3 rmats dataframes and adds their psi values together into 1 dataframe.
    """
    both = trim(df1, ext_label='all events')
    pos_nr = trim(df2, ext_label='nonredundant pos')
    neg_nr = trim(df3, ext_label='nonredundant neg')
    return pd.concat([both, pos_nr, neg_nr])


def violin_plot(df, ax):
    """
    plots violin plot.
    """
    sns.violinplot(x='subset', y=0, data=df, hue='label', split=True, cut=0, ax=ax)


def get_prefix(f):
    """
    returns the prefix of a *.JunctionCountOnly.txt file.
    """
    return os.path.basename(f).split('.')[0]


def make_df_for_violin_plot(f):
    files = get_all_rmats(f)
    if not any(x is None for x in files):
        df = add_labels_and_concat(files[0], files[1], files[2])
        return df
    else:
        return None  # TODO make compatible for < 3 lists


def make_violin_plot_single_rbp(f, output_file):
    """
    from a single original *.JunctionCountOnly.txt file, plot a violinplot of
    density distributions of psi scores for all reps/conditions of that file.
    """
    df = make_df_for_violin_plot(f)
    fig, ax = plt.subplots()
    violin_plot(df, ax=ax)
    plt.legend()
    plt.suptitle(get_prefix(f) + " (p < 0.05, fdr < 0.1, sep > 0.05)")
    plt.tight_layout()
    plt.savefig(output_file)


def get_hist(df, col='IncLevel2', bins=25):
    df['avg'] = df.apply(get_avg_inclevels, axis=1, args=(col,))
    n, bins, patches = plt.hist(df['avg'], bins=bins)
    return n


def make_violin_plots_for_all_jxc_in_directory(input_directory, output_directory, ext='.png'):
    all_files = glob.glob(os.path.join(input_directory, '*.JunctionCountOnly.txt'))

    progress = trange(len(all_files))
    for f in all_files:
        output_file = os.path.join(output_directory, os.path.basename(f).replace('.txt', ext))
        make_violin_plot_single_rbp(f, output_file)
        progress.update(1)


def make_avg_inclevel_dataframe_from_directory_of_jxc(directory, n_threshold=100):
    all_pos_files = glob.glob(os.path.join(directory, '*.JunctionCountOnly.positive.nr.txt'))
    x = defaultdict()
    progress = trange(len(all_pos_files))
    for f in all_pos_files:
        df = pd.read_table(f)
        if df.shape[0] >= n_threshold:
            x[get_prefix(f)] = pd.Series(get_hist(df))
        progress.update(1)
    return pd.DataFrame(x).T


def run_make_inclevel_dataframes():
    pass


def run_num_differential_events(
        clip_manifest_df, rnaseq_manifests_dict, annotation_dir,
        pos_suffix, neg_suffix
):
    pos = defaultdict(dict)
    neg = defaultdict(dict)
    tot = defaultdict(dict)

    print("Missing from folder:\t{}\t{}\t{}\t{}".format(
        "uID", "CELL", "RBP", "REASON")
    )
    for uid in clip_manifest_df['uID']:
        rbpname, cell, _, _, _, rnaseq_expt_prefix, \
        annotation_file = m.get_all_info_for_one_expt(
            uid, annotation_dir, clip_manifest_df, rnaseq_manifests_dict
        )
        if annotation_file is None:
            # Is it because we are setting cutoffs too strict
            # or that we don't have samples downloaded yet?
            print("Missing from folder:\t{}\t{}\t{}\t{}".format(
                uid, cell, rbpname, rnaseq_expt_prefix)
            )
        else:
            df = pd.read_table(annotation_file)
            tot["{}-{}-{}".format(uid, rbpname, cell)] = [
                rnaseq_expt_prefix, df.shape[0]
            ]
            positive, negative = m.get_annotations_from_splicing_prefix(
                annotation_dir,
                rnaseq_expt_prefix,
                pos_suffix,
                neg_suffix
            )
            if positive is None:
                # Not enough significant positive events
                pos["{}-{}-{}".format(uid, rbpname, cell)] = [
                    rnaseq_expt_prefix, pdf.shape[0]
                ]
                print("Missing from folder:\t{}\t{}\t{}\t{}".format(
                    uid, cell, rbpname, 'NO POS SIG EVENTS')
                )
            else:
                pdf = pd.read_table(positive)
                positive_clean = positive.replace(
                    '.MATS.JunctionCountOnly.positive.nr.txt', ''
                )
                pos["{}-{}-{}".format(uid, rbpname, cell)] = [
                    os.path.basename(positive_clean), pdf.shape[0]
                ]
            if negative is None:
                # Not enough significant negative events
                neg["{}-{}-{}".format(uid, rbpname, cell)] = [
                    rnaseq_expt_prefix, pdf.shape[0]
                ]
                print("Missing from folder:\t{}\t{}\t{}\t{}".format(
                    uid, cell, rbpname, 'NO NEG SIG EVENTS')
                )
            else:
                pdf = pd.read_table(negative)
                negative_clean = negative.replace(
                    '.MATS.JunctionCountOnly.negative.nr.txt', ''
                )
                neg["{}-{}-{}".format(uid, rbpname, cell)] = [
                    os.path.basename(negative_clean), pdf.shape[0]
                ]
    pos_df = pd.DataFrame(pos).T
    neg_df = pd.DataFrame(neg).T
    tot_df = pd.DataFrame(tot).T

    merged = pd.merge(
        pos_df, neg_df, how='outer', left_index=True, right_index=True
    )
    merged = pd.merge(
        tot_df, merged, how='outer', left_index=True, right_index=True
    )
    merged.columns = [
        'rnaseq_expt', 'total_events',
        'rnaseq_expt1', 'significant_positive_events',
        'rnaseq_expt2', 'significant_negative_events'
    ]
    del merged['rnaseq_expt1']
    del merged['rnaseq_expt2']

    return merged



def main(argv=None):  # IGNORE:C0111

    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_version_message = '%%(prog)s %s' % (
        program_version,
    )

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument(
        "--clip-manifest",
        dest="clip_manifest",
        required=True,
        help='eclip input normed manifest',
    )
    parser.add_argument(
        "--hepg2-rnaseq-manifest",
        dest="hepg2_rnaseq_manifest",
        required=True,
        help='xintao from graveley lab-styled rnaseq manifest',
    )
    parser.add_argument(
        "--k562-rnaseq-manifest",
        dest="k562_rnaseq_manifest",
        required=True,
        help='xintao from graveley lab-styled rnaseq manifest',
    )
    parser.add_argument(
        "--annotation-parent-dir",
        dest="annotation_dir",
        required=True,
        help='directory where annotations are kept'
    )
    parser.add_argument(
        "--annotation-sub-dir",
        dest="annotation_sub_dir",
        required=True,
        help='subdirectory (i named them as /se/ /a3ss/ etc.'
    )
    parser.add_argument(
        "--pos-suffix",
        dest="pos_suffix",
        required=True,
        help='default: ".positive.nr.txt'
    )
    parser.add_argument(
        "--neg-suffix",
        dest="neg_suffix",
        required=True,
        help='default: ".negative.nr.txt'
    )
    parser.add_argument(
        "--out-file",
        dest="out_file",
        required=True,
        help='out file'
    )
    # Process arguments
    args = parser.parse_args()

    clip_manifest = args.clip_manifest
    hepg2_rnaseq_manifest = args.hepg2_rnaseq_manifest
    k562_rnaseq_manifest = args.k562_rnaseq_manifest
    annotation_dir = args.annotation_dir
    event = args.annotation_sub_dir
    pos_suffix = args.pos_suffix
    neg_suffix = args.neg_suffix
    out_file = args.out_file

    # set variables
    clip_manifest = pd.read_table(clip_manifest)
    rnaseq_manifests = {
        'HepG2': hepg2_rnaseq_manifest,
        'K562': k562_rnaseq_manifest
    }
    annotation_dir = os.path.join(annotation_dir, event)

    # run program
    df = run_num_differential_events(
        clip_manifest, rnaseq_manifests, annotation_dir,
        pos_suffix, neg_suffix
    )

    # write to file
    df.to_csv(out_file, sep='\t')

if __name__ == '__main__':
    main()
