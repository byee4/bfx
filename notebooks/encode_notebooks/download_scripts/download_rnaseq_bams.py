'''
Created on May 19, 2016

@author: brianyee
'''

import sys
import os
import pandas as pd
import urllib, json, requests


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

def main(argv=None): # IGNORE:C0111
    
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    """
    infile: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    Also see: cookbook_notebooks/read_20160408_ENCODE_MASTER_ID_LIST_AllDatasets.ipynb
    """
    parser.add_argument("-i", "--infile", dest="infile", help="input manifest file", required = True )
    parser.add_argument("-o", "--outfile", dest="outfile", help="outfile download list file", required = True )
    parser.add_argument("-l", "--label", dest="label", help="lab name", required = False, default = 'encode-processing-pipeline')
    parser.add_argument("-c", "--cell", dest="cell", help="cell line", required = True, default = 'HepG2')
    
    
    host = 'https://www.encodeproject.org'
    experiments = "https://www.encodeproject.org/experiments/"
    
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    labname = args.label
    cell = args.cell
    
    X = pd.read_table(infile,index_col=0)
    existing = list()
    links = list()
    controls = dict()
    for i in X[X['Cell_line']==cell]['BAM'].dropna():
    # for i in X[X['CellLine']==cell]['ENCODE_ID or GEO Accession'].dropna():  # column value changed when???
    # for i in X[X['CellLine']==cell]['RNASEQ_ENCODEAccID'].dropna():
        print(i)
        url = experiments+i+"/?format=json"
        response = urllib.urlopen(url)
        data = json.loads(response.read())
        if 'code' in data.keys():
            next
        else:
            existing.append(i)
            for i in range(0,len(data['files'])):
                if ((host+data['files'][i]['href']).endswith('bam')) & (data['files'][i]['lab'][u'name'] == labname):
                    #if str(host+data['files'][i]['href']).split('/')[6] in existing_bams:
                    #    continue
                    #else:
                    links.append(host+data['files'][i]['href'])
                    controls[data['accession']] = data['possible_controls'][0]['accession']
    f = open(outfile,'w')
    for i in links:
        f.write("%s\n" % i)
    f.close()
if __name__ == '__main__':
    main()
