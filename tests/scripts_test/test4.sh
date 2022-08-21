#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 4 : Without periodic conditions. Give bad results except if the cluster is far from the box edges."
clustsize   -i data/ -m "H01" -e 100 --mins 5 --eps 3 \
            -o outputs/test4 -vv -f "md.xtc"