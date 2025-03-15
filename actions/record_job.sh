#!/bin/bash

# Get current directory
current_dir=$(pwd)

# Get current date (YYYY-MM-DD)
current_date=$(date '+%Y-%m-%d')

# Append the current directory and date to ~/job_list.txt
echo "$current_date : $current_dir" >> ~/job_list.txt

# Optional: Print a message to indicate success
echo "Recorded $current_dir on $current_date to ~/job_list.txt"
