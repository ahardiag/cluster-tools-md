#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 1 : Use a subset of atoms to apply clustering, e.g. excluding hydrogen atoms"
clustsize   -i data/ -s "resname H01 and not name H* " \
            -e 100 --mins 5 --eps 7.2 \
            -o outputs/test1 -vv -f "md.xtc" \
            --with-double-labelling \
            --ignore-warnings