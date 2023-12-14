#!/usr/bin/python

# Author: Christian Reisswig 2008
# email:  reisswig@aei.mpg.de


import sys
import os
import glob
from numpy import *

usage_msg = """
Extract modes from a decomp_vars output file created by HarmonicDecomposition or SphericalHarmonics.

Usage: sort_hd.py <var_name> <var_nbr> <lmax> <sphere_nbr> <yes|no> <infile> <dir1> ... <dirN>

where:
    <var_name> is a string which will be used to name the output files
    <var_nbr> is the integer 'vars' parameter for the desired variable
    <lmax> is the maximum l-mode which has been extracted
    <sphere_nbr> is the integer number of the extraction sphere
    <yes|no> decides whether to treat this as a complex variable or not (basically we assume that var_nbr+1 is the imaginary part)
    <infile> the filename
    <dirN> are directories containing the decomp_vars* data

Example:
  sort_hd.py psi4 0 4 0 yes harmonicdecomposition::decomposed_vars.xy.asc checkpoint0 checkpoint1 checkpoint2
"""


is_complex = 0


if len(sys.argv)<8 or sys.argv[1] in ['-h', '--help']:
    sys.exit(usage_msg)
var_name = sys.argv[1]
var_nbr, lmax, sphere_nbr = map(int, sys.argv[2:5])
dir_list = sys.argv[7:]
infile = sys.argv[6]

if (sys.argv[5] == "yes"): is_complex = 1

col_ix = 6 #2  # 6
col_iy = 7 #3  # 7
col_time = 9 #4 # 9
col_data = 13 #5 # 13

psi4 = []
psi4_2 = []

last_time = []

print "Processing..."

file_list = []
for d in dir_list:
    file_list.append(str(d)+'/'+infile)  #+thornname+'::decomposed_vars.xy.asc')
    #psi4.set_modes_from_files(file_list, sphere_nbr, var_nbr)
    

for f in file_list:
	sys.stdout.write("Reading "+str(f)+"\n")
	if (os.access(str(f), os.F_OK) == 0):
	    sys.stdout.write("    File does not exist! skipping!\n")
	    continue
        fp = open(str(f), 'r')
        for line in fp:
	    if len(line.split()) == 0 or line[0] == '#':
	        continue 
	    data = []
	    data = map(float, line.split())
	    # make sure we only consider the correct sphere:
	    if (data[col_iy-1] != sphere_nbr):
	        continue
		
	    if len(psi4) <= int(data[col_ix-1]):
		psi4.append([])
		psi4_2.append([])
		last_time.append(-10.0)	
	    time = data[col_time-1]
	    if (last_time[int(data[col_ix-1])] < time):    
		psi4[int(data[col_ix-1])].append([time, complex(data[col_data-1+2*var_nbr], data[col_data-1+2*var_nbr+1])])
		if is_complex:
		    psi4_2[int(data[col_ix-1])].append([time, complex(data[col_data-1+2*(var_nbr+1)], data[col_data-1+2*(var_nbr+1)+1])])
	        last_time[int(data[col_ix-1])] = time	
	fp.close()


print "Number of modes = ", len(psi4), " -> lmax = ", int(sqrt(len(psi4))-1)
print "T_start = ", psi4[0][0][0]
print "T_end = ", psi4[0][-1][0]

if lmax != int(sqrt(len(psi4))-1):
    lmax = int(sqrt(len(psi4))-1)


lm = 0
for l in range(0, lmax+1):
    for m in range(-l, l+1):
	fout = open(var_name+"_l"+str(l)+"m"+str(m)+".dat", 'w')
	#fout = open(var_name+"_"+str(l*l+m)+".dat", 'w')
	for ii in range(0, len(psi4[lm])):
	    if is_complex:
		fout.write(str(psi4[lm][ii][0])+" "+str((psi4[lm][ii][1]+complex(0.0,1.0)*psi4_2[lm][ii][1]).real)+" "+str((psi4[lm][ii][1]+complex(0.0,1.0)*psi4_2[lm][ii][1]).imag)+"\n")
	    else:
		fout.write(str(psi4[lm][ii][0])+" "+str(psi4[lm][ii][1].real)+" "+str(psi4[lm][ii][1].imag)+"\n")
	fout.close()
	lm = lm + 1

print 'Done.'

