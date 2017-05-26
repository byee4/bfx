#!/usr/local/bin/python2.7
# encoding: utf-8
'''
python3.make_23mer_index -- shortdesc

python3.make_23mer_index is a description

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
import re

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2015-08-10'
__updated__ = '2015-08-10'

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
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('-p', "--path", dest="path", help="path to files")
        parser.add_argument('-o', "--output", dest="output", help="output tabbed file")
        parser.add_argument('-xm', "--filterxm", action='store_true', help="filter out XM_ (predicted) genes")
        parser.add_argument('-xr', "--filterxr", action='store_true', help="filter out XR_ (predicted) genes")
        parser.add_argument('-xp', "--filterxp", action='store_true', help="filter out XP_ (predicted) genes")

        # Process arguments
        args = parser.parse_args()

        paths = args.path
        tabbed = args.output
        
        fil = []
        
        if(args.filterxm == True):
            fil.append('ref|XM_')
        if(args.filterxr == True):
            fil.append('ref|XR_')
        if(args.filterxp == True):
            fil.append('ref|XP_')
    
        if(os.path.isfile(tabbed)):
            print("Warning! Output file exists and will be overwritten")
            f = open(tabbed,'w')
            f.close()
        
        print("path to refseq mRNA files: {0}".format(paths))
        seq = ""
        number_of_files = 0
        for name in glob.glob(paths+'*.fna'):
            number_of_files = number_of_files + 1
            print("Adding file: {0}".format(name))
            with open(name, 'r') as infile, open(tabbed, 'a') as outfile:
                head = infile.readline()
                # head = re.findall("ref\|(.*)\|",head)[0].replace('\n','')
                for line in infile:
                    if line[0] != '>':
                        seq = seq + line.replace('\n','')
                    else:
                        for i in range(0,len(seq)):
                            subseq = seq[i:i+23]
                            if subseq[0] == 'G' and subseq[21:23] == "GG":
                                print(subseq[0])
                                if not any(prefix in head for prefix in fil):
                                    outfile.write(">{0}\n{1}\n".format(head,subseq))
                        seq = ""
                        head = line.replace('\n','')
                        head = re.findall("ref\|(.*)\|",head)[0]
                    
        return 0
        
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())