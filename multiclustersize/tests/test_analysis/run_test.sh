#!/bin/bash

# Use a subset of atoms to apply clustering, e.g. excluding hydrogen atoms
../../run_Clustsize.py    -i ../../../traj_tests/ -s "resname H01 and not name H* " \
                    -e 100 --mins 5 --eps 7.2 \
                    -o outputs -vv -f "md.xtc" \
                    --with-double-labelling \
                    --ignore-warnings

exit
# Method 1 : Full Distance Matrix
../../run_Clustsize.py    -i ../../../traj_tests/ -m "H01" -e 100 --mins 5 --eps 3 \
                    -o outputs -vv -f "md.xtc" \
                    --with-distance-matrix

# Method 2 : Double Labelling
../../run_Clustsize.py    -i ../../../traj_tests/ -m "H01" -e 100 --mins 5 --eps 3 \
                    -o outputs -vv -f "md.xtc" \
                    --with-double-labelling \
                    --ignore-warnings

# Method 3 : No periodic conditions
../../run_Clustsize.py    -i ../../../traj_tests/ -m "H01" -e 100 --mins 5 --eps 3 \
                    -o outputs -vv -f "md.xtc"
