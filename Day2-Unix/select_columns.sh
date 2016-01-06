#!/bin/bash

l=1
while IFS=',' read -a line || [[ -n "$line" ]]; do
  if [ $l -eq 1 ]
  then
    perl -le 'print join ",",@ARGV' "${line[1]}" "${line[0]}" "${line[2]}" "${line[3]}"
  else
    if [ ${line[0]} -ge 90 ]
    then
      perl -le 'print join ",",@ARGV' "${line[1]}" "${line[0]}" "${line[2]}" "${line[3]}"
      #echo "${line[@]}"
    fi
  fi
  (( l++ ))
done < "$1"