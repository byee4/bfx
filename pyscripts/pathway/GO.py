#!/usr/bin/env python

"""
GO analysis script
__author__ = 'mlovci'
"""
from __future__ import division
from itertools import izip
from functools import partial

import numpy as np
from scipy.stats import hypergeom
import pandas as pd
import seaborn as sns
import argparse


def hypergeometric(row, lenAllGenes, lenTheseGenes):
    if (row['inBoth'] <= 3) or (row['expressedGOGenes'] < 5):
        return np.nan
    else:
        return hypergeom.sf(row['inBoth'], lenAllGenes,
                            row['expressedGOGenes'], lenTheseGenes)


class GO(object):
    def __init__(self, GOFile):

        # self.GO_to_ENSG = pd.read_table(GOFile, compression="gzip").dropna()
        self.GO_to_ENSG = pd.read_table(GOFile).dropna()

        GO, allGenes = self._generateOntology()
        self.GO = GO
        self.allGenes = allGenes
        self.gene_id_to_name = dict(izip(self.GO_to_ENSG['Gene stable ID'],
                                         self.GO_to_ENSG[
                                             'Gene name']))

    def enrichment(self, geneList, background=None):
        if background is None:
            background = self.allGenes
        return self.GO_enrichment(geneList, expressedGenes=background)

    def _generateOntology(self):
        """
        :return: returns dict of ontologeis and all genes in all go ontologeis as a background
        """
        allGenesInOntologies = set(self.GO_to_ENSG['Gene stable ID'])
        ontology = self.GO_to_ENSG.groupby("GO term accession")
        ontology = ontology.aggregate(lambda x: set(x))
        ontology['nGenes'] = ontology['Gene stable ID'].apply(len)

        return ontology, allGenesInOntologies

    def GO_enrichment(self, geneList, expressedGenes=None):
        geneList = set(list(geneList))
        expressedGenes = set(list(expressedGenes))

        lenAllGenes = len(expressedGenes)
        lenTheseGenes = len(geneList)

        df = self.GO.copy()

        if lenTheseGenes > lenAllGenes:
            raise ValueError(
                "Length of genes examined should not be larger than the total number of genes in organism")

        df['inBoth'] = df['Gene stable ID'].apply(lambda x: len(x & geneList))
        df['expressedGOGenes'] = df['Gene stable ID'].apply(
            lambda x: len(x & expressedGenes))

        hypergeometric_partial = partial(hypergeometric,
                                         lenAllGenes=lenAllGenes,
                                         lenTheseGenes=lenTheseGenes)

        df['Hypergeometric p-Value'] = df.apply(hypergeometric_partial, axis=1)
        num_tests = len(df['Hypergeometric p-Value'].dropna())
        df['Bonferroni-corrected Hypergeometric p-Value'] = df[
            'Hypergeometric p-Value'].apply(
            lambda x: min(x * num_tests, 1)).dropna()

        # Compute various value for backwards compatabality
        df['Ensembl Gene IDs in List'] = df['Gene stable ID'].apply(
            lambda x: x & geneList)
        df['Gene symbols in List'] = df['Ensembl Gene IDs in List'].apply(
            lambda x: {self.gene_id_to_name[gene_id] for gene_id in x})

        df['Ensembl Gene IDs in List'] = df['Ensembl Gene IDs in List'].apply(
            ",".join)
        df['Gene symbols in List'] = df['Gene symbols in List'].apply(",".join)
        df['GO Term Description'] = df['GO term name'].apply(",".join)
        df['GO domain'] = df['GO domain'].apply(",".join)
        df['N Genes in GO category'] = df['Gene stable ID'].apply(len)

        # Rename stuff for backwards compatabality
        df['N Expressed Genes in GO Category'] = df['expressedGOGenes']
        df['N Genes in List and GO Category'] = df['inBoth']

        # Sort
        df.sort_values('Bonferroni-corrected Hypergeometric p-Value', inplace=True)
        # Reorder for presentation purposes
        df = df[[
            'GO Term Description',
            'Bonferroni-corrected Hypergeometric p-Value',
            'N Genes in List and GO Category',
            'N Expressed Genes in GO Category',
            'N Genes in GO category',
            'Ensembl Gene IDs in List',
            'Gene symbols in List',
            'GO domain',
            'Hypergeometric p-Value',
        ]]
        return df

    @staticmethod
    def enrichment_score_vectorized(hit_values, miss_values):
        normalized_hit_values = hit_values / hit_values.sum()
        normalized_miss_values = miss_values / miss_values.sum()

        enrichment_score = normalized_hit_values.cumsum() - normalized_miss_values.cumsum()
        return enrichment_score

    @staticmethod
    def get_largest_score_vectorized(enrichment_score):
        combined_scores = pd.concat({"max_score": enrichment_score.max(),
                                     "min_score": enrichment_score.min()}).unstack()
        return combined_scores.apply(
            lambda x: x.max_score if np.abs(x.max_score) > np.abs(
                x.min_score) else x.min_score, axis=0)

    @staticmethod
    def normalized_enrichment_score_vectorized(normalized_hit_values,
                                               normalized_miss_values):
        enrichment_score = normalized_hit_values.cumsum() - normalized_miss_values.cumsum()
        return enrichment_score

    @staticmethod
    def fast_normalized_enrichment_score_vectorized(normalized_hit_values,
                                                    normalized_miss_values):
        enrichment_score = normalized_hit_values.as_matrix().cumsum(
            axis=0) - normalized_miss_values.as_matrix().cumsum(axis=0)
        return enrichment_score

    @staticmethod
    def fast_get_largest_score_vectorized(enrichment_score):
        # Enrichment score is already a matrix
        # I intentionally ignore WHERE the max score is here for speed purposes.  Its not important for the shuffling, only the non-shuffled case
        max_scores = enrichment_score.max(axis=0)
        min_scores = enrichment_score.min(axis=0)
        max_scores[np.abs(max_scores) < np.abs(min_scores)] = min_scores[
            np.abs(max_scores) < np.abs(min_scores)]
        return max_scores

    def gsea(self, gene_list, max_size=500, min_size=25, num_iterations=1000):

        """
        :param gene_list: pandas series where index is ensembl gene ids and values are scores
        :param max_size: max size of go terms or gene lists to allow into GSEA analysis
        :param min_size: min size of go terms of gene lists to allow into GSEA analysis
        :param num_iterations: number of random iterations to perform
        :return:
        """

        gene_list = gene_list.sort_values(ascending=False)
        self.gene_list = gene_list
        # get proper sets
        large_sets = self.GO[(self.GO['nGenes'] > min_size) & (
        self.GO['nGenes'] < max_size)].copy()

        # set up matix of genes in each set
        result = {}
        in_set = {}
        for name, genes in large_sets['Gene stable ID'].iteritems():
            result[name] = gene_list.copy()
            in_set[name] = {gene_id: True for gene_id in genes}

        # fill out set values with the rest of the genes in the gene list
        result = pd.DataFrame(result)
        in_set = pd.DataFrame(in_set)
        in_set = in_set.fillna(False)

        incomplete_in_set = in_set.T

        for gene in result.index.difference(in_set.index):
            incomplete_in_set[gene] = False

        in_set = incomplete_in_set[result.index].T

        # create matrix of genes in each go term and not in each go term (hit and miss)
        hit_values = result.copy()
        miss_values = result.copy()

        hit_values[~in_set] = 0
        hit_values = np.abs(hit_values)
        self.hit_values = hit_values

        miss_values[:] = 0
        miss_values[~in_set] = 1

        # calculate enrichment scores for true values and store for plotting later
        enrichment_score = self.enrichment_score_vectorized(hit_values,
                                                            miss_values)
        self.enrichment_score = enrichment_score

        largest_enrichment = self.get_largest_score_vectorized(
            enrichment_score)

        # normalize hit and miss values first for speed
        normalized_hit_values = hit_values / hit_values.sum()
        normalized_miss_values = miss_values / miss_values.sum()

        # generate random list
        results = {}
        for x in range(num_iterations):
            index = list(hit_values.index)
            np.random.shuffle(index)
            shuffled_hit_values = normalized_hit_values.ix[index]
            shuffled_miss_values = normalized_miss_values.ix[index]

            results[x] = self.fast_get_largest_score_vectorized(
                self.fast_normalized_enrichment_score_vectorized(
                    shuffled_hit_values, shuffled_miss_values))

        shuffled_results = pd.DataFrame(results)
        shuffled_results.index = hit_values.columns
        shuffled_results = shuffled_results.T

        # Compute p-values as in paper, generated z-scores based on only positive or negative distributions

        pos_enrichment = largest_enrichment[largest_enrichment >= 0]
        neg_enrichment = largest_enrichment[largest_enrichment < 0]

        pos_shuffled = shuffled_results[pos_enrichment.index]
        neg_shuffled = shuffled_results[neg_enrichment.index]

        pos_shuffled = pos_shuffled[pos_shuffled >= 0]
        neg_shuffled = neg_shuffled[neg_shuffled < 0]

        pos_count = pos_shuffled.count()
        neg_count = neg_shuffled.count()

        pos_p_value = (pos_shuffled > pos_enrichment).sum() / pos_count
        neg_p_value = (neg_shuffled < neg_enrichment).sum() / neg_count
        p_values = np.abs(pd.concat([pos_p_value, neg_p_value]))

        pos_nes_enrichment = pos_enrichment / pos_shuffled.mean()
        neg_nes_enrichment = (neg_enrichment / neg_shuffled.mean()) * -1

        pos_nes_shuffled = pos_shuffled / pos_shuffled.mean()
        neg_nes_shuffled = (neg_shuffled / neg_shuffled.mean()) * -1

        pos_p_value_nes = (
                          pos_nes_shuffled > pos_nes_enrichment).sum() / pos_count
        neg_p_value_nes = (
                          neg_nes_shuffled < neg_nes_enrichment).sum() / neg_count

        nes_p_values = np.abs(pd.concat([pos_p_value_nes, neg_p_value_nes]))
        nes_enrichment = pd.concat([pos_nes_enrichment, neg_nes_enrichment])
        # format output

        enrichment_df = pd.concat({"enrichment": largest_enrichment,
                                   'nes_enrichment': nes_enrichment,
                                   "p-Value": p_values,
                                   'nes_p-Value': nes_p_values
                                   }).unstack().T

        enrichment_df = enrichment_df.join(large_sets)

        num_tests = len(p_values.dropna())
        enrichment_df['Bonferroni-corrected p-Value'] = enrichment_df[
            'p-Value'].apply(lambda x: min(x * num_tests, 1)).dropna()
        enrichment_df['Bonferroni-corrected NES p-Value'] = enrichment_df[
            'nes_p-Value'].apply(lambda x: min(x * num_tests, 1)).dropna()

        enrichment_df['GO Term Description'] = enrichment_df[
            'GO Term Name'].apply(",".join)
        enrichment_df['GO domain'] = enrichment_df['GO domain'].apply(",".join)

        enrichment_df = enrichment_df[[
            'GO Term Description',
            'enrichment',
            'nes_enrichment',
            'nes_p-Value',
            'Bonferroni-corrected NES p-Value',
            'p-Value',
            'Bonferroni-corrected p-Value',
            'Ensembl Gene ID',
            'Gene name',
            'Ensembl Transcript ID',
            'GO Term Name',
            'GO Term Definition',
            'GO Term Evidence Code',
            'GO domain',
            'GOSlim GOA Description',
            'GOSlim GOA Accession(s)',
            'nGenes'
        ]]

        return enrichment_df, hit_values, enrichment_score


def plot_go_term(go_term, gene_list, hit_values, enrichment_score, fig):
    gene_list = gene_list.sort_values(ascending=False)

    genes_in_set = pd.DataFrame(hit_values[go_term].copy())
    genes_in_set['x_loc'] = xrange(len(genes_in_set))
    genes_in_set['y_loc'] = genes_in_set[go_term].apply(
        lambda x: 1 if x > 0 else 0)
    genes_in_set = genes_in_set[genes_in_set.y_loc == 1]

    ax = fig.add_subplot(3, 1, 1)
    ax.scatter(genes_in_set.x_loc, genes_in_set.y_loc, s=5, alpha=.6)
    ax.set_xlim(0, len(enrichment_score[go_term]))
    sns.despine(ax=ax, left=False, bottom=True)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_ylabel("Gene in set")

    ax = fig.add_subplot(3, 1, 2)
    ax.plot(gene_list)
    ax.set_xlim(0, len(enrichment_score[go_term]))
    sns.despine(ax=ax, bottom=True)
    ax.set_xticklabels([])
    ax.set_ylabel("Correlation")

    ax = fig.add_subplot(3, 1, 3)
    ax.plot(enrichment_score[go_term])
    sns.despine(ax=ax)
    ax.set_xlim(0, len(enrichment_score[go_term]))
    ax.set_ylabel("Enrichment Score")
    ax.set_xlabel("Gene Rank")


def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--GO_file",
        required=True,
    )
    parser.add_argument(
        "--genes_list",
        required=True,
    )
    parser.add_argument(
        "--bg_list",
        required=True,
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    args = parser.parse_args()
    go_file = args.GO_file
    genes_list = args.genes_list
    bg_list = args.bg_list
    out_file = args.out_file

    go = GO(go_file)
    df = go.enrichment(genes_list, background=bg_list)
    df.to_csv(
        out_file
    )
if __name__ == "__main__":
    main()
