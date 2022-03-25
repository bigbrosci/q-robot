#!/usr/bin/env bash 
# To save the CONTCAR for specific ionic steps 
if [ -e CONTCAR ]; then 
    step=$(grep F= OSZICAR | tail -n 1 | awk '{print $1}')
    #echo $step
    cp CONTCAR CONTCAR-$step
fi 
