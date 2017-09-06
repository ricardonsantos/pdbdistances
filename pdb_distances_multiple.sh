#!/bin/bash
FILES=*.pdb
for i in $FILES
do
   ./pdb_calpha_dist.exe $i
done
