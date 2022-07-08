#!/usr/bin/env bash 
SLEEP=60
path=/home/2092/job_infor
for i in {1..10000} ;  do
echo 'Check' "$i" 
python3 job_check.py
sleep $SLEEP
done

