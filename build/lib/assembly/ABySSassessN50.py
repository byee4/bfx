#!/usr/local/bin/python2.7
# encoding: utf-8
'''
assembly.soap_avg_coverage -- shortdesc

assembly.soap_avg_coverage is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os
import glob

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2015-04-27'
__updated__ = '2015-04-27'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2015 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-i", "--input", help="directory")
        max_contig = 0
        max_scaffold = 0
        longest_contig = 0
        longest_scaffold = 0
        # Process arguments
        args = parser.parse_args()
        for file in glob.glob(args.input+"/*/*"):
            if file.endswith(".csv"):
                f = open(file,'r')
                lines = f.readlines()
                for line in lines:
                    line = line.split(',')
                    print(file.replace(args.input,'')+":\t"+line[5]+"\t"+line[8]+"\t"+line[10]),
                    if (str.isdigit(line[5])):
                        n50 = int(line[5])  # n50 value
                        long = int(line[8]) # max contig/scaff/unitig size
                        if ('contigs' in line[10]):
                            if(n50 > max_contig):
                                max_contig = n50
                                max_contig_file = file.replace(args.input,'')
                            if( long > longest_contig):
                                longest_contig = long
                                longest_contig_file = file.replace(args.input,'')
                        if ('scaffolds' in line[10]):
                            if(n50 > max_scaffold):
                                max_scaffold = n50
                                max_scaffold_file = file.replace(args.input,'')
                            if( long > longest_scaffold):
                                longest_scaffold = long
                                longest_scaffold_file = file.replace(args.input,'')
        print("------------------------------------------------------------------------")
        print("highest N50 value (contig): {0} found in: {1}".format(max_contig, max_contig_file))
        print("highest N50 value (scaffold): {0} found in: {1}".format(max_scaffold, max_scaffold_file))
        print("longest contig size: {0} found in: {1}".format(longest_contig, longest_contig_file))
        print("longest scaffold size: {0} found in: {1}".format(longest_scaffold, longest_scaffold_file))
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())