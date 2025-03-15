#!/usr/bin/env bash 
tar -zcf "$1".tar "$1" && rm "$1" -fr 
