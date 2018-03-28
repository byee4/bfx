__author__ = 'brian'
import matplotlib
matplotlib.use('Agg')
import glob
import pandas as pd
import seaborn as sns
import os
from argparse import ArgumentParser
from tqdm import trange
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import rmats_inclevel_analysis as rmats


def get_merged_avg_jc(jc_list, min_events=0):
    merged = pd.DataFrame()
    progress = trange(len(jc_list))
    for f in jc_list:
        df = pd.read_table(f)
        if df.shape[0] >= min_events:
            name = os.path.basename(f).replace('SE.MATS.JunctionCountOnly','').replace('.txt','')
            df['junction'] = rmats.get_junction_region(df)
            dx = rmats.get_avg_dpsi_for_all_junctions(df)
            dx.columns = [name]
            merged = pd.merge(merged, dx, how='outer', left_index=True, right_index=True)
        progress.update(1)
    return merged


def plot_correlation(jc_list, out_file, method='spearman'):
    merged = get_merged_avg_jc(jc_list)
    g = sns.clustermap(merged.fillna(0).corr(method=method), xticklabels=False, figsize=(10,30))
    x = plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)  # For y axis
    plt.savefig(out_file)
    return merged


def main():
    parser = ArgumentParser()

    parser.add_argument(
        "-o", "--output",
        dest="output",
        help="output svg",
        required=True
    )
    parser.add_argument(
        "-i", "--input",
        dest="input",
        help="input folder",
        required=True
    )
    parser.add_argument(
        "-k", "--keep",
        action='store_true',
        help="keep intermediates",
        default=False,
        required=False
    )
    parser.add_argument(
        "-s", "--suffix",
        dest="suffix",
        help="default: \'SE.MATS.JunctionCountOnly.txt\'",
        default='-SE.MATS.JunctionCountOnly.txt',
        required=False
    )
    parser.add_argument(
        "-m", "--method",
        dest="method",
        help="correlation: [spearman], pearson, kendall",
        default='spearman',
        required=False
    )

    args = parser.parse_args()

    i = args.input
    o = args.output
    s = args.suffix
    k = args.keep
    m = args.method

    all_files = glob.glob(os.path.join(i,'*{}'.format(s)))
    print(len(all_files))
    df = plot_correlation(all_files, o, m)

    if k:
        df.to_csv(os.path.splitext(o)[0] + '.txt', sep='\t')

if __name__ == "__main__":
    main()
