#!/bin/bash

echo -e "\n----------------------------> TEST <--------------------------------"
echo "test 5 : plot results of test 1"
plotclustersize -i outputs/test1 -o outputs/test1/figures/ --tmax 100 --rename
