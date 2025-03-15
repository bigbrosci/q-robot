#!/bin/bash

# Extract and display only the job paths from ~/job_list.txt
awk -F' : ' '{print $2}' ~/job_list.txt

# Optional: Print a message if no jobs are found
if [ ! -s ~/job_list.txt ]; then
  echo "No jobs recorded in ~/job_list.txt"
fi
