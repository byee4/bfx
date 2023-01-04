__author__ = 'gpratt'

import gzip
from optparse import OptionParser
import sys

barcodes = """AAGCAAT A01
GGCTTGT B06
ACAAGTT C01
TGGTCCT D08fixed
ATGACCNNNNT  A03
TCCTGTNNNNT  G07
CAGCTTNNNNT  A04
GGATACNNNNT  F05
NNNNNCCTATAT X1A
NNNNNTGCTATT X1B
NNNNNTATACTT X2A
NNNNNATCTTCT X2B""".split("\n")


barcodes = dict([item.split() for item in barcodes])
barcode_name_to_sequence = {value.strip(): key.strip() for key, value in barcodes.items()}

def add_back_barcodes(in_file_name, out_file_name, barcode_name):
    #reads through initial file parses everything out
    with gzip.open(in_file_name) as fastq_file, gzip.open(out_file_name, 'w') as out_file:
        while True:
            try:
                name_1 = fastq_file.next()
                seq_1 = fastq_file.next()
                fastq_file.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality_1 = fastq_file.next()

                randomer = name_1.split(":")[0][1:]
                name_1 = "@" + ":".join(name_1.split(":")[1:])
                seq_1 = barcode_name_to_sequence[barcode_name] + seq_1
                quality_1 = len(barcode_name_to_sequence[barcode_name]) * "-" + quality_1

                out_file.write(name_1)
                out_file.write(seq_1)
                out_file.write(plus)
                out_file.write(quality_1)

            except StopIteration:
                break

def add_back_randomers(in_file_name, out_file_name):
    #reads through initial file parses everything out
    with gzip.open(in_file_name) as fastq_file, gzip.open(out_file_name, 'w') as out_file:
        while True:
            try:
                name_1 = fastq_file.next()
                seq_1 = fastq_file.next()
                fastq_file.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality_1 = fastq_file.next()

                randomer = name_1.split(":")[0][1:]
                name_1 = "@" + ":".join(name_1.split(":")[1:])
                seq_1 = randomer + seq_1
                quality_1 = len(randomer) * "-" + quality_1

                out_file.write(name_1)
                out_file.write(seq_1)
                out_file.write(plus)
                out_file.write(quality_1)

            except StopIteration:
                break

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--in_file", help="bam file to add barcodes back to")
    parser.add_option("-o", "--out_file", help="out bam file")
    parser.add_option("-b", "--barcode", help="barcode (if used)", default=None)
    (options, args) = parser.parse_args()

    if options.barcode:
        add_back_barcodes(options.in_file, options.out_file, options.barcode)
    else:
        add_back_randomers(options.in_file, options.out_file)

    sys.exit(0)
