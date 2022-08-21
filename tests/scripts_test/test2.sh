#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 2 : Full Distance Matrix (SLOW but exact)"
clustsize   -i data/ -m "H01" -e 100 --mins 5 --eps 3 \
            -o outputs/test2 -vv -f "md.xtc" \
            --with-distance-matrix