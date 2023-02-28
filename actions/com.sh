#!/usr/bin/env bash
file=$1
find  -name "$1" -execdir  tar  -zcvf   "$1".tar   --remove-files {}   \+
