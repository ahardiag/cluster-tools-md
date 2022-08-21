#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 3 : Double Labelling with all atoms and a cutoff of 3A"
clustsize   -i data/ -m "H01" -e 100 --mins 5 --eps 3 \
            -o outputs/test3 -vv -f "md.xtc" \
            --with-distance-matrix