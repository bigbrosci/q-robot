#!/usr/bin/env bash 
for file in LOCPOT_*; do if [[ ! $file == *.tar ]]; then tar -zcvf "$file.tar" "$file" && rm -f "$file"; fi; done
