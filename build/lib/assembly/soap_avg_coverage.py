#!/usr/local/bin/python2.7
# encoding: utf-8
'''
python3.soap_avg_coverage -- shortdesc

python3.soap_avg_coverage is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os
import re

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
        parser.add_argument("-i", "--input", help=".contig file")

        # Process arguments
        args = parser.parse_args()
        i = 0
        sum_cov = 0
        bases = 0
        with open(args.input, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    regexp = re.compile('cvg_(\d*\.{1}\d*)')
                    cov = regexp.findall(line)
                    sum_cov = sum_cov + float(cov[0])
                    i = i + 1
                else:
                    line = line.replace("\n","")
                    bases = bases + len(line)
        print("avg coverage = {0:.2f}".format(sum_cov/i))
        print("avg depth = {0:.2f}".format(bases))
        print("3,532,909")
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())