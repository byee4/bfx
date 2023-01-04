#!/bin/bash

awk '$4>=3 && $5>=3' ${1}
