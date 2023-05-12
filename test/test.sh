#!/bin/bash

dir1="out"
dir2="output_sequences"

differing_files=()
num_differing_files=0

for i in {63..123}
do
    file1="${dir1}/node_$i.fas"
    file2="${dir2}/node_$i.fas"

  if [ -f "$file1" ] && [ -f "$file2" ]; then
    if ! diff -q "$file1" "$file2" > /dev/null; then
      differing_files+=("$i.fas")
      ((num_differing_files++))
    fi
  fi
done

echo "Number of differing files: $num_differing_files"
echo "Differing files:"
for file in "${differing_files[@]}"; do
  echo "$file"
done