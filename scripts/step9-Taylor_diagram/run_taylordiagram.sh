#!/bin/bash

module load ncl

ncl 1-find_statistics.ncl >log
sed -e s///g log >log2
sort log2>log3

rm -f log log2

# 1. You will need to vi log3, move the title line to the top and remove the other unwanted lines
#    after that, use command:
#              sed -e s/\(0\)\	//g log3 >statistics_mass_Tbio.txt
#
# 2. Then rearrange and extract data from statistics_mass_Tbio.txt
#    ratio.txt (from column 5), corr.txt  (from column 6)
#    Both text files have the format:
#                        ---  Model -->
#      |  Forest 
#      |  fraction
#      v
#
# 3. You have the inputs for Taylor diagrams now!
#
# 4. Original NCL script of Taylor diagram can be found at NCAR NCL website
#    https://www.ncl.ucar.edu/Applications/taylor.shtml
