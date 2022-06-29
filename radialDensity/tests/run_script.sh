#!/bin/bash

# Compute the radial density of all three components of the system
# around the 10-molecules aggregates
../radial_density_aggregates.py  -i ../../traj_tests/ \
                                -c ../../multiclustersize/tests/test_analysis/outputs/ \
                                -r init.gro -f md.xtc \
                                -s Data_H01_100_wat-EtOH_m5_box10 \
                                -z 10
