#!/usr/env python

"""
See original notebook for interactive usage:
pilot_notebooks/pilot_create_non_overlapping_rmats_junction_intervals.ipynb
"""

import subset_rmats_junctioncountonly as s

import pybedtools

def fake_completely_overlapping_bedtool_as_df():
    intervals = []
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '500', '2500', 'big', '0.75', '+'
            ]))
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '1000', '2000', 'small', '0.5', '+'
            ]))
    df = pybedtools.BedTool(intervals).to_dataframe()
    df.columns = ['chr','exonStart_0base','exonEnd','geneSymbol','IncLevelDifference','strand']
    return df

def fake_partially_overlapping_bedtool_as_df():
    intervals = []
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '500', '1000', 'genex', '0.75', '+'
            ]))
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '750', '1500', 'genex', '0.5', '+'
            ]))
    df = pybedtools.BedTool(intervals).to_dataframe()
    df.columns = ['chr','exonStart_0base','exonEnd','geneSymbol','IncLevelDifference','strand']
    return df

def fake_share_start_exon_bedtool_as_df():
    intervals = []
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '500', '1000', 'genex', '0.75', '+'
            ]))
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '500', '1500', 'genex', '0.5', '+'
            ]))
    df = pybedtools.BedTool(intervals).to_dataframe()
    df.columns = ['chr','exonStart_0base','exonEnd','geneSymbol','IncLevelDifference','strand']
    return df

def fake_share_end_exon_bedtool_as_df():
    intervals = []
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '500', '1500', 'genex', '0.75', '+'
            ]))
    intervals.append(pybedtools.create_interval_from_list([
                'chr1', '750', '1500', 'genex', '0.5', '+'
            ]))
    df = pybedtools.BedTool(intervals).to_dataframe()
    df.columns = ['chr','exonStart_0base','exonEnd','geneSymbol','IncLevelDifference','strand']
    return df


def test_create_non_overlapping_regions_from_rmats_df_1():
    print("check if appropriately completely overlapping regions are properly split and scores averaged")
    print("original regions: "
          "[chr1, 500, 2500, 0.75, big, -],"
          "[chr1, 1000, 2000, 0.5, small, -]")
    df = fake_completely_overlapping_bedtool_as_df()

    new_df = s.create_non_overlapping_regions_from_rmats_df(df)

    assert new_df['start'][0] == 500
    assert new_df['end'][0] == 1000
    assert new_df['score'][0] == 0.750

    assert new_df['start'][1] == 1000
    assert new_df['end'][1] == 2000
    assert new_df['score'][1] == 0.625

    assert new_df['start'][2] == 2000
    assert new_df['end'][2] == 2500
    assert new_df['score'][2] == 0.750


def test_create_non_overlapping_regions_from_rmats_df_2():
    print("exon1 start < exon2 start, exon1 end < exon2 end")
    print("original regions: "
          "[chr1, 500, 1000, 0.75, big, +],"
          "[chr1, 750, 1500, 0.5, small, +]")
    df = fake_partially_overlapping_bedtool_as_df()

    new_df = s.create_non_overlapping_regions_from_rmats_df(df)

    assert new_df['start'][0] == 500
    assert new_df['end'][0] == 750
    assert new_df['score'][0] == 0.750

    assert new_df['start'][1] == 750
    assert new_df['end'][1] == 1000
    assert new_df['score'][1] == 0.625

    assert new_df['start'][2] == 1000
    assert new_df['end'][2] == 1500
    assert new_df['score'][2] == 0.500


def test_create_non_overlapping_regions_from_rmats_df_3():
    print("Share starting exon position")
    print("original regions: "
          "[chr1, 500, 1000, 0.75, genex, +],"
          "[chr1, 500, 1500, 0.5, genex, +]")
    df = fake_share_start_exon_bedtool_as_df()

    new_df = s.create_non_overlapping_regions_from_rmats_df(df)

    assert new_df['start'][0] == 500
    assert new_df['end'][0] == 1000
    assert new_df['score'][0] == 0.625

    assert new_df['start'][1] == 1000
    assert new_df['end'][1] == 1500
    assert new_df['score'][1] == 0.500


def test_create_non_overlapping_regions_from_rmats_df_4():
    print("Share ending exon position")
    print("original regions: "
          "[chr1, 500, 1500, 0.75, genex, +],"
          "[chr1, 750, 1500, 0.5, genex, +]")
    df = fake_share_end_exon_bedtool_as_df()

    new_df = s.create_non_overlapping_regions_from_rmats_df(df)

    assert new_df['start'][0] == 500
    assert new_df['end'][0] == 750
    assert new_df['score'][0] == 0.750

    assert new_df['start'][1] == 750
    assert new_df['end'][1] == 1500
    assert new_df['score'][1] == 0.625



