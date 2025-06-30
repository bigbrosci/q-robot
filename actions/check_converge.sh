#!/bin/bash

# File to check
OUTCAR_FILE=$1

# Check if the job finished by searching for "Voluntary"
if tail -n 10 "$OUTCAR_FILE" | grep -q "Voluntary"; then
    JOB_FINISHED="Yes"
else
    JOB_FINISHED="No"
fi

# Extract NSW value
NSW=$(grep "NSW" $OUTCAR_FILE | head -n 1 | awk '{print $3}')

# Extract NELM value
NELM=$(grep "NELM" $OUTCAR_FILE | head -n 1 | awk '{print $3}')

# Extract the last iteration details (ionic and electronic steps)
LAST_ITER=$(grep "Iter" $OUTCAR_FILE  | tail -n 1 | sed 's/-//g' | sed 's/Iteration//g')
IONIC_STEP=$(echo $LAST_ITER | awk -F'(' '{print $1}')
ELECTRONIC_STEP=$(echo $LAST_ITER | awk -F'(' '{print $2}' | tr -d ')')

echo $OUTCAR_FILE
if [[ $JOB_FINISHED == "Yes" ]]; then
    if [[ $NSW -eq 0 || $NSW -eq 1 ]]; then
        if [[ $NELM -gt $ELECTRONIC_STEP ]]; then
            echo "$OUTCAR_FILE, Good: Job finished. Ionic: $IONIC_STEP, NSW: $NSW, Electronic: $ELECTRONIC_STEP, NELM: $NELM." >> check_results.out
            echo "$OUTCAR_FILE" >> list_good.txt
        else
            echo "$OUTCAR_FILE, Bad: Job finished but not converged. Ionic: $IONIC_STEP, NSW: $NSW, Electronic: $ELECTRONIC_STEP, NELM: $NELM." >> check_results.out
            echo "$OUTCAR_FILE" >> list_bad.txt
        fi
    else
        if [[ $NSW -gt $IONIC_STEP && $NELM -gt $ELECTRONIC_STEP ]]; then
            echo "$OUTCAR_FILE, Good: Job finished. Ionic: $IONIC_STEP, NSW: $NSW, Electronic: $ELECTRONIC_STEP, NELM: $NELM." >> check_results.out
            echo "$OUTCAR_FILE" >> list_good.txt
        else
            echo "$OUTCAR_FILE, Bad: Job finished but not converged. Ionic: $IONIC_STEP, NSW: $NSW, Electronic: $ELECTRONIC_STEP, NELM: $NELM." >> check_results.out
            echo "$OUTCAR_FILE" >> list_bad.txt
        fi
    fi
else
    echo "$OUTCAR_FILE, Bad: Job did NOT finish." >> check_results.out
    echo "$OUTCAR_FILE" >> list_bad.txt
fi

