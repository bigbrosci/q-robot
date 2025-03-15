for i in OUTCAR*; do echo $i $(grep '  without' $i|tail -n 1|awk '{print $NF}'); done
