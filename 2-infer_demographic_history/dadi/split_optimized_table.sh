#!/bin/bash

# define model
model=($(echo "${1}"))
# define pops used
set=($(echo "${2}"))

# queries line numbers
rep_lines=($(grep -n "AIC" ${set}.${model}.optimized.txt | cut -d':' -f1))
# array numbers
rep_len=${#rep_lines[@]}

for (( i=1; i<${rep_len}+1; i++ ))
 do  
  # define lines of repetition results
  if [[ ${i} -lt ${rep_len} ]]
   then 
    START=($(echo "${rep_lines[${i}-1]}"))
    END=($(echo "${rep_lines[${i}]} - 1" | bc -l))
  elif [[ ${i} -eq ${rep_len} ]]
   then
    START=($(echo "${rep_lines[${i}-1]}"))
    END=($(cat ${set}.${model}.optimized.txt | wc -l))
  fi
  # extract lines of query results
  echo "rep $i starts at $START and ends at $END"
  cat ${set}.${model}.optimized.txt | sed -n "${START},${END}p" > ${set}.${model}.optimized.r_${i}.txt
done
