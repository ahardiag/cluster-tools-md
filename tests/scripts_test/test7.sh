#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 7 : compute eccentricity for clusters found during test 1"
eccentricity    -i data/ \
                -c outputs/test1 \
                -r init.gro -f md.xtc \
                -s Data_H01_100_wat-EtOH_m5_box10 \
                --outfigdirname outputs/test1/figures