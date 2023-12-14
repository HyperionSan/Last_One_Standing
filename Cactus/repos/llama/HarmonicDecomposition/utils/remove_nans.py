#!/usr/bin/python2.4

# Author: Christian Reisswig 2008
# email: reisswig@aei.mpg.de


import sys
import os
import glob


usage_msg = """
USAGE: remove_nans.py <file>

Removes NaNs from any ASCII-datafile and replaces them by 0.

"""

if (len(sys.argv) < 2):
    sys.exit(usage_msg)    


filename=sys.argv[1]


fp = open(filename, 'r')
fout = open(filename+".no_NANs", 'w')
for line in fp:
    if line[0]=='#' or line[0]=='"' or len(line.strip())==0:
	fout.write(line)
	continue
    data = []
    
    data = line.split()
    for ii in range(0, len(data)):
	if (data[ii].strip() == 'nan') : data[ii] = "0"
	fout.write(data[ii].strip()+" ")
    fout.write("\n")
    
    
    
fp.close()
fout.close()

