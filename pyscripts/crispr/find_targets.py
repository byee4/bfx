'''
Created on Aug 21, 2015

@author: brianyee
'''
import subprocess
import re
import math

from general import seq

'''MAC'''
bowtie = "/Applications/bowtie-1.1.1/bowtie"
idx = "/Users/brianyee/Documents/Datasets/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set"
biomart_file = "/Users/brianyee/Documents/Datasets/mart_export.txt"

'''Linux
bowtie = "/usr/bin/bowtie"
idx = ""
biomart_file = "/home/brian/Documents/CRISPR/mart_export.txt"
'''


def run_bowtie(bowtie_loc, bowtie_idx, sequence):
    """ 
    runs bowtie according to preset parameters according to: 
    http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0124633

    The "core" includes the last 12nt proximal to PAM.
    -e is set to 150 to 'force' bowtie into finding 4+ mismatches.

    Actual parameters are: -a, -n 3, -l 12, -e 150

    args: 
        bowtie: location of bowtie
        idx: index of genome 
        sequence: query sequence

    returns: 
        http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-outputs
        the line-delimited output of bowtie *not including*:
            name of alignment (col 1)
            read quality (col 6)
            number of other instances aligned (col 7)
    """
    cmd = [
        bowtie_loc, bowtie_idx, '--quiet', '-p', '4', '-n', '3',
        '-l', '12', '-e', '150', '--suppress', '1,6,7', '-y', '-a',
           '--chunkmbs', '256', '-c', sequence
    ]
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
    output = output.split('\n')
    return output


def mit_evaluate(sequence):
    """ 
    takes an entire sequence and filters out only those 23mers
    containing NGG, and evaluates each 
    args: 
        sequence: the entire sequence being evaluated

    prints: 
        a list of candidate sequences and their aggregate scores
    """
    count = 0  # number of candidate sequences
    candidates = {}
    for h in range(0, len(sequence) - 23):
        subsequence = sequence[h:h + 23].upper()
        subsequence_rc = seq.revcomp(subsequence)
        if subsequence[21:23] == "GG" and subsequence not in candidates:
            candidates[subsequence] = mit_score(subsequence)
            print("{0} has a score of: {1} with a hit count of {2}".format(
                subsequence, candidates[subsequence][0],
                candidates[subsequence][1])
            )
            count = count + 1
        if subsequence_rc[21:23] == "GG" and subsequence_rc not in candidates:
            candidates[subsequence_rc] = mit_score(subsequence_rc)
            count = count + 1
            print(
            "{0} has a score of: {1} with a hit count of {2}".format(
                subsequence_rc, candidates[subsequence_rc][0],
                candidates[subsequence_rc][1])
            )
    for key in candidates:
        print("{0} -> {1}".format(key, candidates[key]))

    print("{0} candidates evaluated.".format(count))


def mit_score(subsequence):
    """ computes and returns the aggregate score for each candidate sequence 
        The formula is found at: http://crispr.mit.edu/about
        The higher the score, the better the candidate

        args: 
            subsequence: candidate 23mer with PAM sequence

        returns: 
            float value pertaining to score
            number of hits bowtie found (might/might not be important?)
    """
    output = run_bowtie(subsequence)

    hit_count = 0  # number of bowtie hits that are less or equal to 4 MM
    exact_matches = 0  # we are hoping for just one match in genome, the one specified by query

    all_scores = []  # all the scores taken from each hit of bowtie
    for i in range(0, len(output) - 1):

        output[i] = output[i].split('\t')

        mismatch = re.sub("[ATCG:>]", "", output[i][4]).split(',')
        if (not mismatch[0]):  # we have encountered an exact match
            exact_matches = exact_matches + 1
            all_scores.append(1)
        else:  # not counting the exact hits
            pos = map(float, mismatch)  # numerizes mismatch string to float
            if len(pos) <= 4:  # if 4 or less mismatches (sometimes bowtie will return a ton of mismatches)
                all_scores.append(mit_single_hit_score(pos))
                hit_count = hit_count + 1
    aggregate_score = (100.0) / (100.0 + sum(all_scores) - 1)
    return aggregate_score, hit_count


def mit_single_hit_score(mismatch_pos=[], *args):
    """ returns single hit score of MIT's algorithm 
        found at: http://crispr.mit.edu/about
        The higher the score, the more likely this hit 
        will produce an off target effect

        args: 
            mismatch_pos[]: array containing mismatch positions

        returns: 
            float value pertaining to hit score
    """
    score1 = 1.0
    # print(mismatch_pos)
    matrix = [0, 0, 0.014, 0,
              0, 0.395, 0.317, 0,
              0.389, 0.079, 0.445, 0.508,
              0.613, 0.851, 0.732, 0.828,
              0.615, 0.804, 0.685, 0.583]  # pre-computed weights for mismatch positions
    for i in range(0, len(mismatch_pos)):
        score1 = score1 * (1.0 - matrix[i])

    davg = (max(mismatch_pos) - min(mismatch_pos)) / (len(mismatch_pos) - 1)
    score2 = 1.0 / ((19 - davg) / 19) * 4 + 1
    score3 = 1.0 / (len(mismatch_pos) * len(mismatch_pos))
    return score1 * score2 * score3


def cctop_evaluate(sequence):
    """ takes an entire sequence and filters out only those 23mers
        containing NGG, and evaluates each 
        args: 
            sequence: the entire sequence being evaluated

        prints: 
            a list of candidate sequences and their aggregate scores
    """
    count = 0  # number of candidate sequences
    candidates = {}
    for h in range(0, len(sequence) - 23):
        subsequence = sequence[h:h + 23].upper()
        subsequence_rc = seq.revcomp(subsequence)
        if subsequence[21:23] == "GG" and subsequence not in candidates:
            candidates[subsequence] = cctop_score(subsequence)
            # print("{0} has a score of: {1} with a hit count of {2}".format(subsequence, candidates[subsequence][0], candidates[subsequence[1]]))
            count = count + 1
        if subsequence_rc[21:23] == "GG" and subsequence_rc not in candidates:
            candidates[subsequence_rc] = cctop_score(subsequence_rc)
            count = count + 1
            # print("{0} has a score of: {1} with a hit count of {2}".format(subsequence_rc, candidates[subsequence_rc][0], candidates[subsequence_rc[1]]))
    for key in candidates:
        print("{0} -> {1}".format(key, candidates[key]))

    print("{0} candidates evaluated.".format(count))


def cctop_score(subsequence):
    """ computes and returns the aggregate score for each candidate sequence 
        The formula is found at: http://crispr.mit.edu/about
        The higher the score, the better the candidate

        args: 
            subsequence: candidate 23mer with PAM sequence

        returns: 
            float value pertaining to score
            number of hits bowtie found (might/might not be important?)
    """
    output = run_bowtie(subsequence)  # Note: this runs bowtie against the ENTIRE sequence (including the PAM)
    hit_count = 0  # number of bowtie hits that are less or equal to 4 MM
    exact_matches = 0  # we are hoping for just one match in genome, the one specified by query
    total_score = 0  # total score for the subsequence

    offtargets = []  # array of dictionaries of potential offtarget candidates

    for i in range(0, len(output) - 1):  # for each hit

        dist = 100000  # distance from the nearest exon (max 100k)
        output[i] = output[i].split('\t')  # parse output

        sequence = subsequence  # starting value is the query sequence
        strand = output[i][0]  # strand either '+' or '-'
        chromosome = re.sub("chr", "", output[i][1])  # ie "chr1 -> 1"
        start = int(output[i][2])  # start position of hit
        mismatch = output[i][4]  # mismatch "code" provied by bowtie's default output (offset:reference>read,)
        mismatch_positions = re.sub("[ATCG:>]", "", output[i][4]).split(
            ',') if mismatch else []  # array of mismatch positions only

        if (not mismatch):  # we have found an exact hit
            exact_matches = exact_matches + 1
        if (len(mismatch_positions) < 5):  # we are not interested in 5+ mismatched hits
            if (strand == '-'):
                sequence = seq.mutate_sequence(sequence, mismatch)
                # sequence = seq.revcomp(seq.get_sequence_from_das("hg38",chromosome,start-2,start+20)).upper() # let's change this to use bowties mismatches later.
            else:
                sequence = seq.mutate_sequence(sequence, mismatch)
                # sequence = seq.get_sequence_from_das("hg38",chromosome,start+1,start+23).upper()
            if (sequence[21:23] == "GG"):
                hit_count = hit_count + 1
                dist = seq.get_closest_exon(chromosome, start, biomart_file)
                if (0 < dist[0] < 100000):  # if location start is intronic or intergenic and has associated exon
                    logdist = math.log10(dist[0])
                else:
                    logdist = 0  # If target site and exon coordinates overlap, the distance is assigned to 0
                offtargets.append(
                    {"Sequence": sequence, "MM": mismatch_positions, "Distance": logdist, "Gene": dist[1]})
    for j in range(0, len(offtargets)):
        score = (offtargets[j]["Distance"] + cctop_single_hit_score(offtargets[j]["MM"])) / len(offtargets)
        print("offtarget sequence: {0}\tmismatch positions: {1}\tclosest associated gene: {2}\tscore: {3}".format(
            offtargets[j]["Sequence"], offtargets[j]["MM"], offtargets[j]["Gene"], score))
        total_score = total_score + score

    print("candidate has {0} exact matches with {1} offtarget hits with 4 or less mismatches.".format(exact_matches,
                                                                                                      len(offtargets)))
    print("score before off_target subtraction: {0}".format(total_score))
    return total_score - (len(offtargets) - 1)


def cctop_single_hit_score(mismatch_pos=[], *args):
    """ returns single hit score of CCTop's algorithm 
        found at: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0124633
        The higher the score, the more likely this hit 
        will produce an off target effect

        args: 
            mismatch_pos[]: array containing mismatch positions

        returns: 
            float value pertaining to hit score
    """
    score = 0
    if (len(mismatch_pos) == 0):
        return 1

    for i in range(0, len(mismatch_pos)):
        score = score + pow(1.2, int(mismatch_pos[i]))
    return score


def main():
    print(cctop_score("GCAAAGGAAATTGACAGCAATGG"))  # GGGACACTTTGCGTTCGGGCTGG
    # print(g.revcomp(g.get_sequence_from_das("hg38","chr17",7687420-2,7687420+20)))
    # print(g.revcomp(g.get_sequence_from_das("hg38","chr17",72519513-2,72519513+20)))
    # print(g.get_sequence_from_das("hg38","chr1",12007429+1,12007429+23))


if __name__ == '__main__':
    main()