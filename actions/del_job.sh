#!/bin/bash

# Get the current directory
current_dir=$(pwd)

# Remove the current directory from job_list.txt
sed -i "\|$current_dir|d" ~/job_list.txt

# Optional: Print a message indicating success
echo "Removed $current_dir from ~/job_list.txt"

