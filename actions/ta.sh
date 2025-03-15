#!/usr/bin/env bash 
for i in * ; do 
	if [ -e $i/OUTCAR ]; then 
		e=$(grep '  without' $i/OUTCAR |tail -n 1 | awk '{print $7}')
		echo $i' '$e
	fi
done 
