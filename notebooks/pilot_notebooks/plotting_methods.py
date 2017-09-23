# from clip_analysis.src import regions
# from clip_analysis.src import parsers as p
import parsers as p
import matplotlib.pyplot as plt
from itertools import izip
import seaborn as sns
import numpy as np
import os

# do this to avoid making tons of dots in svg files:
import matplotlib
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

palette = sns.color_palette("hls", 8)

REGIONS = {
    'noncoding_exon':palette[0],
    '3utr':palette[1],
    '5utr':palette[2],
    'intron':palette[3],
    'noncoding_intron':palette[4],
    'CDS':palette[5],
    'intergenic':palette[6],
    '5utr_and_3utr':palette[7]
}


def plot_ip_foldchange_over_input_reads(
        ip_l2fc, inp_reads_by_loc,
        out_file,
        field_list=REGIONS,
        alpha=0.3
):
    """
    Plots the region-based analysis of genes enriched in IP over
     reads in size-matched input. This is the same figure used in Figure2b
     of Ann's IMP paper.

    Parameters
    ----------
    ip_l2fc : basestring
    inp_reads_by_loc : basestring
    out_file : basestring
        output file
    field_list : dict
        dictionary of {region: color} to plot
    Returns
    -------

    """
    df = p.scatter_matrix(ip_l2fc, inp_reads_by_loc)

    f, ax = plt.subplots(figsize=(10, 10))
    for region, color in field_list.iteritems():
        if region in df.columns:
            ax.scatter(
                np.log2(df[region] + 1),
                df["{} l2fc".format(region)],
                color=color,
                alpha=alpha
            )
        else:
            print("region {} not found in dataframe matrix".format(
                region
            ))
    ax.set_title("Region-based analysis of genes enriched")
    ax.set_xlabel("Reads in SMInput (log2)")
    ax.set_ylabel("Fold Enrichment (log2)")
    plt.legend()
    plt.savefig(out_file)


def plot_region_distribution(
        wd, out_file, l10p, l2fc,
        l10p_col=3, l2fc_col=4,
        trim_suffix=".peaks.l2inputnormnew.bed.compressed.bed.annotated",
        format="eric"
):
    """
    Takes a directory containing '.annotated' bedfiles and plots a
    stacked barchart of the peaks within regions (defined by REGIONS).
    This will also, as a side effect, produce filtered bedfile for each
    globbed *.annotated

    Parameters
    ----------
    wd : basestring
    out_file : basestring
        output file
    l10p : float
        log10 p value
    l2fc : float
        log2 fold change
    trim_suffix : basename
        suffix to trim from the filename (makes plot easier to read)

    Returns
    -------

    """
    out_dir = os.path.dirname(out_file)
    cumsum_events = p.get_cumulative_sum_counts(
        wd, l10p, l2fc, l10p_col, l2fc_col, out_dir, trim_suffix, format, REGIONS
    )

    f, ax = plt.subplots(figsize=(10, 10))

    legend_builder = []
    legend_labels = []
    for region, color in izip(
            reversed(cumsum_events.index),
            sns.color_palette("Set2", len(cumsum_events.index))
    ):
        if trim_suffix != None:
            cumsum_events.columns = [
                event.replace(trim_suffix, '') for event in
                cumsum_events.columns
                ]

        names = np.array(
            ["".join(item) for item in cumsum_events.columns]
        )

        sns.barplot(
            names,
            y=cumsum_events.ix[region], color=color, ax=ax
        )

        legend_builder.append(
            plt.Rectangle((0, 0), .25, .25, fc=color, edgecolor='none')
        )
        legend_labels.append(region)

    sns.despine(ax=ax, left=True)

    ax.set_ylim(0, 1)

    l = ax.legend(legend_builder,
                  legend_labels, loc=1, ncol=1,
                  prop={'size': 12},
                  bbox_to_anchor=(1.4, 0.8))
    l.draw_frame(False)
    [tick.set_rotation(90) for tick in ax.get_xticklabels()]

    ax.set_ylabel("Fraction of Peaks", fontsize=14)
    [tick.set_fontsize(12) for tick in ax.get_xticklabels()]
    ax.set_title(
        "Fraction of Peaks among RBPs \
        \n-log10(p-value):{}, log2(fold-change):{}".format(
            l10p, l2fc)
    )
    f.savefig(out_file)