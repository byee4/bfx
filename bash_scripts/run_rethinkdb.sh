#!/bin/bash

source activate rethinkdb;

RETHINKDB_PORT=8579

rethinkdb --http-port $RETHINKDB_PORT

ssh -NR $RETHINKDB_PORT:localhost:$RETHINKDB_PORT tscc-login1 &
