#!/bin/bash

for file in cice6_*CICE6-*.png; do
    newfile="$file"
    newfile="${newfile/Jan/01}"
    newfile="${newfile/Feb/02}"
    newfile="${newfile/Mar/03}"
    newfile="${newfile/Apr/04}"
    newfile="${newfile/May/05}"
    newfile="${newfile/Jun/06}"
    newfile="${newfile/Jul/07}"
    newfile="${newfile/Aug/08}"
    newfile="${newfile/Sep/09}"
    newfile="${newfile/Oct/10}"
    newfile="${newfile/Nov/11}"
    newfile="${newfile/Dec/12}"
    mv "$file" "$newfile"
done
