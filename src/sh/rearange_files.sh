#!/bin/bash

for file in cice6_??_*.png; do
    month="${file:6:2}"
    region="${file##*_}"   # Extracts the region (assuming it's the last part before .png)
    region="${region%.png}"  # Remove the .png extension
    mkdir -p "$region/$month"
    #echo "moving $file to $region/$month/"
    mv "$file" "$region/$month/"
done

