#!/usr/bin/env python

"""
This script will compare two bed files to see how similar they are.
From two bed files, this script will first merge overlapping regions,
then intersects the two bed files. Intersecting peaks that overlap by
at least 50% of each other are counted, and the number of intersecting
regions divided by the total number of regions for each bed file
will provide a similarity score
"""

# transitionning to python2/python3 support
# uncomment from this compatibility import list, as py3/py2 support progresses
from __future__ import print_function
from __future__ import division
# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
# from future.utils import iteritems

import glob
import os
import argparse
import pprint
import pymongo


def mongodb_search(user, password, nickname, summary):
    client = pymongo.MongoClient(
        "mongodb://{}:{}@cluster0-shard-00-00-hcum8.mongodb.net:27017,cluster0-shard-00-01-hcum8.mongodb.net:27017,cluster0-shard-00-02-hcum8.mongodb.net:27017/test?ssl=true&replicaSet=Cluster0-shard-0&authSource=admin&retryWrites=true".format(
            user,
            password
        )
    )
    db = client['u19-test']
    # db['genomics'].create_index([('summary', 'text')])
    if nickname is not None:
        for x in (db['genomics'].find({"aggr_nickname": nickname})):
            pprint.pprint(x)
    if summary is not None:
        for x in (db['genomics'].find({"$text": {"$search": summary}},{'_txtscr': {'$meta': 'textScore'}}).sort([('_txtscr', {'$meta':'textScore'})])):
            pprint.pprint("################################################################################")
            pprint.pprint(x)


def main():
    parser = argparse.ArgumentParser(description="Compares two bed files")
    parser.add_argument("--password", help="password", required=True, default="")
    parser.add_argument("--user", help="user", required=True, default="bay001")
    parser.add_argument("--nickname", help="search experiment", required=False, default=None)
    parser.add_argument("--summary", help="search summary", required=False, default=None)
    args = parser.parse_args()

    mongodb_search(args.user, args.password, args.nickname, args.summary)


if __name__ == '__main__':
    main()