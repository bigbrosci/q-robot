#!/usr/bin/env bash 

for i in $(cat list); do                                                                    
if [ -e $i/OUTCAR ]; then                                                                   
  echo  -e $i "\t" $(grep '  without' $i/OUTCAR |tail -n 1 |awk '{print $7}')               
else                                                                                        
  echo -e $i "\t"                                                                           
fi                                                                                          
done                                                                                        
