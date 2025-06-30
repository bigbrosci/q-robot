#!/bin/bash

# enter_job.sh — Enter the working directory of a Slurm job

if [ "$#" -ne 1 ]; then
    echo "Usage: ej.sh <job_ID>"
    exit 1
fi

job_info=$(scontrol show job "$1")
work_dir=$(echo "$job_info" | grep -oP 'WorkDir=\K\S+')

if [ -d "$work_dir" ]; then
    cd "$work_dir" || exit
    echo "✅ Entered directory: $work_dir"
    # Optional: start interactive shell
    exec $SHELL
else
    echo "❌ Directory not found: $work_dir"
    exit 2
fi

