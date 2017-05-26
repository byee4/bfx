#!/usr/local/bin/python2.7
# encoding: utf-8
'''
rnaseq.mapping -- shortdesc

rnaseq.mapping is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os

from general import functions
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from subprocess import Popen
from subprocess import PIPE

__all__ = []
__version__ = 0.1
__date__ = '2015-09-17'
__updated__ = '2015-09-17'

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
        #parser.add_argument('-t', '--tophat', action='store_true', dest='tophat', default=False)
        parser.add_argument('-t', action='store_true', default=True, dest='tophat', help='Use the tophat program')
        parser.add_argument('-i', action='store', dest='index', help='path to Bowtie index')
        # Process arguments
        args = parser.parse_args()
        bowtie_index = args.index
        
        if(args.tophat):
            tophat_path = functions.which("tophat")
            cmd = [tophat_path, bowtie_index]
            print(cmd)
            output = Popen(cmd, stdout=PIPE).communicate()[0]

        
        return 0
    
    
        
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    main()
    sys.exit(0)