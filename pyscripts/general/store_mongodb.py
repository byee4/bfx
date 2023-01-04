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
import json

def read_jsonlike(jsonlike_file):
    """
    Reads a json-like (json object with a shebang header) and returns the json object.

    :param jsonlike_file: basestring
    :param logger:
    :return data:
    """
    with open(jsonlike_file, 'r') as f:
        f.readline()  # skips the /usr/bin/env line
        try:
            data = json.load(f)
            # data = json_util.loads(data)
        except ValueError:
            return 1
    data['mdsum'] = os.path.splitext(os.path.basename(jsonlike_file))[0].split('_')[-1]
    print(data['mdsum'])
    return data

def mongodb_insert(user, password, jsonlike_file):
    j = read_jsonlike(jsonlike_file=jsonlike_file)
    client = pymongo.MongoClient(
        "mongodb://{}:{}@cluster0-shard-00-00-hcum8.mongodb.net:27017,cluster0-shard-00-01-hcum8.mongodb.net:27017,cluster0-shard-00-02-hcum8.mongodb.net:27017/test?ssl=true&replicaSet=Cluster0-shard-0&authSource=admin&retryWrites=true".format(
            user,
            password
        )
    )
    db = client['u19-test']
    # db['genomics'].insert_one(j)
    if len(db['genomics'].find_one({"mdsum": j['mdsum']})) == 0:
        db['genomics'].insert_one(j)

def main():
    parser = argparse.ArgumentParser(description="Compares two bed files")
    parser.add_argument("--password", help="password", required=True, default="")
    parser.add_argument("--user", help="user", required=True, default="bay001")
    parser.add_argument("--jsonlike", help="search experiment", required=False, default=None, nargs='+')

    args = parser.parse_args()
    user = args.user
    password = args.password

    for jsonlike_file in args.jsonlike:
        mongodb_insert(
            user=user,
            password=password,
            jsonlike_file=jsonlike_file
        )


if __name__ == '__main__':
    main()