#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 8 : Compute the radial density of all three components of the system
around the 10-molecules aggregates found in test 1."
radialdens  -i data/ \
            -c outputs/test1/ \
            -r init.gro -f md.xtc \
            -s Data_H01_100_wat-EtOH_m5_box10 \
            -z 10 \
            --outfigdirname outputs/test1/figures
