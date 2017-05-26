'''
Created on Apr 24, 2015

@author: http://news.open-bio.org/news/author/peterc/
'''

#This Python script requires Biopython 1.51 or later
from Bio import SeqIO

def interleave(iter1, iter2) :
    for (forward, reverse) in zip(iter1,iter2):
        assert forward.id == reverse.id
        forward.id += "/1"
        reverse.id += "/2"
        yield forward
        yield reverse

def main():
    #Setup variables (could parse command line args instead)
    file_f = "/Users/brianyee/Documents/Datasets/Genomics/run_1.trim.paired.fastq"
    file_r = "/Users/brianyee/Documents/Datasets/Genomics/run_2.trim.paired.fastq"
    file_out = "/Users/brianyee/Documents/Datasets/Genomics/run.interleaved.fastq"
    format = "fastq" #or "fastq-illumina", or "fasta", or ...
    records_f = SeqIO.parse(open(file_f,"rU"), format)
    records_r = SeqIO.parse(open(file_r,"rU"), format)
    
    handle = open(file_out, "w")
    count = SeqIO.write(interleave(records_f, records_r), handle, format)
    handle.close()
    print("%i records written to %s" % (count, file_out))
    
if __name__ == '__main__':
    main()