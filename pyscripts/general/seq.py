"""
These are old scripts that I haven't tested or used in a very long while.
"""
import pymysql
import re
import os
import urllib
import sqlite3
import glob
from Bio import SeqIO

# from urllib.request import urlopen  # Python3 
import xml.etree.ElementTree as etree
from collections import defaultdict


translation_dict_dna={'A':'T','C':'G','T':'A','G':'C','N':'N'}
valid_fasta_file_extensions=['.fa','.fasta']
valid_gb_file_extensions=['.gb','.genbank','.gbk']


def comp(seq,dictionary=translation_dict_dna):
    """ Standard complement sequence function.
    
        Returns: the complement sequence
    """
    seq=seq.upper()
    comp=[None]*len(seq)
    for i in range(0,len(seq)):
        comp[i]=dictionary[seq[i]]
    return ''.join(comp)


def revcomp(seq,dictionary=translation_dict_dna):
    """ Standard reverse complement function.
    
        Returns: the reverse complement
    """
    seq=seq.upper()
    comp=[None]*len(seq)
    for i in range(0,len(seq)):
        comp[len(seq)-1-i]=dictionary[seq[i]]
    return ''.join(comp)


def gc(seq):
    """ Standard GC Content calculator.
    
        Returns: the GC Content (as a percentage integer) of the sequence
    """
    count=0
    for letter in seq:
        if letter.upper()=='C' or letter.upper()=='G':
            count=count+1
    return (count/float(len(seq)))*100


def get_ext(my_file):
    """
    Just helps determine whether or not a file extension is a
    "FASTA" one or a "GENBANK" one.
    Returns: either "fasta" or "genbank"
    """
    if os.path.splitext(my_file)[1] in valid_fasta_file_extensions:
        return 'fasta'
    elif os.path.splitext(my_file)[1] in valid_gb_file_extensions:
        return 'genbank'


def query_genbank(infile):
    """
    Returns a record given a genbank file.
    The function assumes the file is valid and that it has the fields:
        LOCUS, ACCESSION, ORGANISM, CDS, ORIGIN and SEQUENCE
    If these fields are blank or nonexistent, it will return the record
    with a blank for those fields.

    Args:
        infile: The genbank file

    Returns:
        A dictionary (record) containing:
        {locus, accession, organism, features, sequence}
    """
    sequence=""
    locus_id=""
    accession=""
    organism=""
    features=defaultdict(list)
    records=[]
    seq_flag=False
    warn=False
    regex=re.compile("[atcg]",re.IGNORECASE)
    with open(infile) as f:
        for line in f:
            if line[:5]=="LOCUS":
                locus_id=line[12:12+line[12:].find(' ')]
            elif line[:9]=="ACCESSION":
                accession=line[12:12+line[12:].find('\n')]
            elif line[2:10]=="ORGANISM":
                organism=line[12:12+line[12:].find('\n')]
            elif line[5:8]=="CDS":
                cds=[int(line[21:21+line[21:].find('..')]),
                int(line[21+line[21:].find('..')+2:line.find('\n')])]
                features['cds'].append(cds)
            elif line[0:6]=="ORIGIN":
                seq_flag=True
            if seq_flag==True and line[0:6]!="ORIGIN":
                sequence=sequence+''.join(regex.findall(line))
            if line[:2]=="//":
                if accession not in records:
                    records.append({'locus':locus_id,
                                        'accession':accession,
                                        'organism':organism,
                                        'features':features,
                                        'sequence':sequence.upper()})
                else: 
                    warn=True
                features=defaultdict(list)
                sequence=""
                seq_flag=False
    if warn==True:
        print("Warning: duplicate accession numbers found.")
    return records


def get_sequence_from_das(my_db,chromosome,start,end):
    """
    Queries the UCSC DAS server given the database, chromosome and position
    Returns the sequence if it is given valid coordinates, else an empty
    string.

    Args:
        my_db: the database name (can be hg38, etc.)
        chromosome: given as "chr1," "chr2," etc.
        start: the raw start position from the beginning of the chromosome
        end: The raw end position

    Returns:
        Sequence from the DAS

    See:
        http://stackoverflow.com/questions/22328319/python-fetch-sequence-from-das-by-coordinates
    """
    sequence=''
    base = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment="%my_db
    url = base + chromosome + ':' + str(start) + ',' + str(end)
    # xml = urlopen(url).read() # Python3
    xml = urllib.urlopen(url).read() # Python2 module
    if xml != '':
        w = open('temp.xml', 'w')
        w.write(xml)
        w.close()
        doc = etree.parse('temp.xml',parser=etree.XMLParser())
        if doc != '':
            # should only return one sequence, right???
            sequence = doc.findall('./SEQUENCE/DNA')[0].text.replace('\n','')
        else:
            print('THE SEQUENCE DOES NOT EXIST FOR GIVEN COORDINATES')
    return sequence


def query_refseq_ucsc(refseq_query,exon_specific=False,my_db="hg38",my_user="genome",
                      my_host="genome-mysql.cse.ucsc.edu"):
    """
    Returns a record of the gene associated with the refseq accession number
    (ie. NM_000000). The record contains information about the gene, including:
        name, accession, cds_start, cds_end, sequence, strand, chromosome

    Args:
        refseq_query: a refseq accession number string
        exon_specific: sometimes we don't need the sequence, but rather
            the exons to avoid targeting sequences that don't exist prior
            to processing. Setting this flag to True ensures a return of
            an array rather than a string in the dictionary.
        my_db: Defaults to the latest human hg38 database. Can be any
            database in UCSC
        my_user: Defaults to "genome" so unless you're a superstar,
            don't change it.
        my_host: Defaults to "genome-mysql.cse.ucsc.edu."
    Returns:
        A dictionary containing:
            name, accession, cds_start, cds_end, sequence, strand, chromosome
            If "exon_specific" is True, the 'sequence' will be an array
            containing individual exons instead of a sequence string.
    """
    myDB = pymysql.connect(host=my_host,user=my_user,db=my_db)
    cur = myDB.cursor()
    record = dict()
    seq = ""
    refseq_query = refseq_query[:refseq_query.find('.')] # we are unconcerned with versions
    try:
        cur.execute("select r.chrom, r.exonCount, r.exonStarts, r.exonEnds, \
                     r.strand, \
                     c.name, \
                     f.geneName \
                      from cds c\
                       join gbCdnaInfo i on c.id = i.cds \
                       join refGene r on r.name = i.acc \
                       join refFlat f on f.name = i.acc \
                    where i.acc = '%s'"%refseq_query)
        result = cur.fetchone()
        if result is None:
            return -1
        if result[5] != 'n/a':
            cds_start = int(result[5][:result[5].find('..')])
            cds_end = int(result[5][(result[5].find('..')+2):])
        else: 
            cds_start = 0
            cds_end = 0
        num_exons = result[1]
        exon_arr = [0 for i in range(num_exons)]
        exon_start = result[2].split(',')
        exon_end = result[3].split(',')
        strand = 0
        if result[4] == '-':
            strand = -1
            for i in range(0,num_exons):
                exon = revcomp(get_sequence_from_das(my_db,result[0],
                               (int(exon_start[num_exons-(i+1)])+1),
                                exon_end[num_exons-(i+1)]).upper())
                exon_arr[i] = exon
                seq=seq+exon
        else:
            strand = 1
            for i in range(0,num_exons):
                exon = get_sequence_from_das(my_db,result[0],(int(exon_start[i])+1),
                                  exon_end[i]).upper()
                exon_arr[i] = exon
                seq=seq + exon
        if (exon_specific):
            record={'name':result[6],'accession':refseq_query,'cds_start':cds_start,
                'cds_end':cds_end,'chromosome':result[0],
                'sequence':exon_arr,'strand':strand}
        else:
            record={'name':result[6],'accession':refseq_query,'cds_start':cds_start,
                'cds_end':cds_end,'chromosome':result[0],
                'sequence':seq,'strand':strand}
    except IndexError:
        print("cannot find sequence!")
    return record


def get_closest_exon(chromosome,pos1,pos2,db="ensembl.sql"):
    """
    Gets the closest exons (max 100kb away) from the position specified.
    It uses (and creates if necessary)
    """
    
    gene_start = pos1 if pos1 < pos2 else pos2
    gene_end = pos2 if pos2 > pos1 else pos1

    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute("SELECT * FROM exons WHERE chromosome = '{0}' AND \
                ({1} BETWEEN exon_start-100000 AND exon_end+100000 OR \
                {2} BETWEEN exon_start-100000 AND exon_end+100000)".format(
        chromosome,gene_start,gene_end
    ))
    ret = cur.fetchall()
    
    min_dist = 100000
    min_gene = ""
    
    
    for i in range(0,len(ret)):
        exon_start = int(ret[i][4])
        exon_end = int(ret[i][5])
        gene = ret[i][2]
        if((exon_start <= gene_start <= exon_end) or (exon_start <= gene_end <= exon_end)):
            return 0, gene # not SUPER accurate, we avoid alternative transcripts
        else:
            dist = min(abs(gene_end-exon_start),abs(gene_start-exon_end))
            if(dist < min_dist):
                min_dist = dist
                min_gene = gene
    return min_dist, min_gene


def create_refseq_exon_sqlite(infile, db="ensembl.sql"):
    """
    Creates the exon sqlite table from results downloaded from biomart.
    To get this data, go to: http://www.biomart.org/biomart/martview/
    and select your species and dataset. I selected attributes:
        Ensembl Gene ID,
        Exon Chr Start (bp)
        Exon Chr End (bp)
        Associated Gene Name
        Chromosome Name
        Strand

    Args:
        infile: file downloaded (tab-delimited) from biomart that includes
        exon position information
        db: name of the sql file to be created
    """
    con = sqlite3.connect(db)
    cur = con.cursor()
    try:
        # create the exon table
        cur.execute("""CREATE TABLE exons(key INTEGER PRIMARY KEY AUTOINCREMENT, \
                gene_id CHARACTER(20), \
                gene_name TEXT COLLATE NOCASE, \
                chromosome CHARACTER(5), \
                exon_start INTEGER(20), \
                exon_end INTEGER(20), \
                strand INTEGER(1))""")
        con.commit()
        
        with open(infile,'rb') as p:
            for line in p:
                line = line.decode(encoding='UTF-8',errors='strict')
                line = line.replace("\n","").split('\t')
                
                try:
                    cur.execute("""INSERT INTO exons( \
                            gene_id,gene_name,chromosome,exon_start,exon_end,strand) \
                            VALUES(?,?,?,?,?,?)""", \
                            (line[0],line[4],line[3],line[1],line[2],line[5]))
                except sqlite3.IntegrityError:
                    print("choking on line:")
                    print(line)
        con.commit()
        print("Done creating the exon sqlite tables.")
    except sqlite3.OperationalError: 
        pass
        # print("Tables already made!")


def mutate_sequence(sequence,mismatch):
    """
    Returns the reference given mismatches (as formatted by standard
    Bowtie Output). The Bowtie mismatch string follows the
    "offset:reference-base>read-base" format, and assumes the given
    sequence is the read string.

    Args:
        sequence: The unadulterated read sequence
        mismatch: The mismatch pattern

    Returns:
        sequence after mismatch pattern is applied
    """
    if(not mismatch): # if there is no mismatch, return original sequence!
        return sequence

    mismatch = mismatch.split(',') # parse mismatch pattern
    new_sequence = ['']*len(sequence) # empty placeholder for new sequence
    positions = {} # dictionary that ultimately holds offset:character
    
    # find mismatches and add to dictionary
    for i in range(0,len(mismatch)):
        positions[int(mismatch[i][:mismatch[i].find(':')])] = \
        mismatch[i][mismatch[i].find(':')+1:mismatch[i].find('>')]
    
    
    for i in range(0,len(sequence)):
        if i in positions:
            new_sequence[i] = positions[i]
        else:
            new_sequence[i] = sequence[i]
    new_sequence = "".join(new_sequence)
    return new_sequence


def fix_fasta_headers(infile, outfile, regex_pattern):
    """
    Formats a fasta file header using regex pattern.
    So instead of the normal header, new headers are just
    the first found instance of the regex capture.
    
    Parameters
    ----------
    infile : basestring
    outfile : basestring
    regex_pattern : basestring

    Returns
    -------
    
    """
    
    o = open(outfile, 'w')
    with open(infile,'r') as f:
        for line in f:
            if line[0] == '>':
                o.write(">{0}\n".format(re.findall(regex_pattern, line)[0]))
            else:
                o.write(line)


def append_fasta(some_path, out_file, ext='fa'):
    """
    Appends/concatenates fasta files in a path to each other.
    
    Parameters
    ----------
    some_path : basestring
    out_file : basestring
    ext : basestring
        default: 'fa'
    Returns
    -------

    """
    o = open(out_file, 'w')
    for fasta_file in glob.glob(some_path+"*.{}".format(ext)):
        print(fasta_file)
        if fasta_file != out_file:
            with open(fasta_file, 'r') as f:
                for line in f:
                    o.write(line)


def get_seq_dict_from_file(f, seq_ids=[], file_type='fasta', equal=True):
    """
    Returns dictionary of {name : sequence}

    Parameters
    ----------

    f : basestring
        file location of a fasta file
    seq_ids : list
        list of sequence ids to search. Empty field returns all sequences
    equal : bool
        True if seq_id needs to be identical
        False if we just have partial identifier
    file_type : basestring
        default "fasta" type file
    Returns
    -------
    records : dict
        {name:sequence}
    """
    records = {}
    for record in SeqIO.parse(f, file_type):
        if len(seq_ids) > 0:
            for name in seq_ids:
                if equal:
                    if name == record.id:
                        records[record.id] = record.seq
                else:
                    if name in record.id:
                        records[record.id] = record.seq
        else:
            records[record.id] = record.seq
    return records


def get_seq_sizes(f, seq_ids=[], file_type='fasta', equal=True):
    """
    Returns dictionary of {name : seqsize}

    Parameters
    ----------
    f
    seq_ids
    equal
    file_type

    Returns
    -------

    """
    lengths = {}
    records = get_seq_dict_from_file(f, seq_ids, file_type, equal)

    for seq_id, sequence in records.iteritems():
        lengths[seq_id] = len(sequence)
    return lengths


def write_chrom_sizes(
        infile, outfile, seq_ids=[], file_type='fasta', equal=True
):
    """
    Writes a chrom sizes file, useful for any software that requires
    a chrom.sizes tabbed file on custom genomes.
    
    Parameters
    ----------
    infile : basestring
    outfile : basestring
    seq_ids : list
    file_type : basestring
    equal : bool

    Returns
    -------

    """
    chrom_sizes = get_seq_sizes(infile, seq_ids, file_type, equal)
    with open(outfile, 'w') as o:
        for chrom, length in chrom_sizes.iteritems():
            o.write("{}\t{}\n".format(chrom, length))

