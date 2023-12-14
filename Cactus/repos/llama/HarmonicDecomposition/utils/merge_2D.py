#!/usr/bin/python2.4

# Author: Christian Reisswig 2008
# email: reisswig@aei.mpg.de


import sys
import os
import glob
#from numpy import *


usage_msg = """
USAGE: merge.py <directories> <file>

Merges 2D-Carpet ASCII-files.

"""

if (len(sys.argv) < 3):
    sys.exit(usage_msg)    


numberofdirectories = len(sys.argv)-1

curdir = os.path.abspath(os.getcwd())

alldirs = sys.argv[1:numberofdirectories]

os.chdir(alldirs[0])
allfiles = glob.glob(sys.argv[numberofdirectories])
os.chdir("..")
alldata=[]

mergeddir = curdir+"/"+os.path.basename(alldirs[0])+".MERGED"
print "Putting merged data to ", mergeddir

os.mkdir(mergeddir)

os.chdir(alldirs[0]+"/..")

for file in allfiles:
	data = []
	for dir in alldirs:
		data.append([])
		filename=str(dir)+"/"+str(file)
		fp = open(filename, 'r')
		for line in fp:
			#if line[0]=='#' or line[0]=='"' or len(line.strip())==0:
			#	continue
			data[-1].append(line)
		fp.close()
	filename=str(mergeddir+"/"+os.path.basename(file))
	fp = open(filename, 'w')

	for jj in range(0, len(data[0])):
	    #for kk in range(0, len(data[0][jj])):
		#fp.write(data[0][jj][kk]+" ")
	    #fp.write("\n")
	    fp.write(data[0][jj])
	    if (len(data[0][jj].split()) > 0 and data[0][jj][0] != '#'):
		lasttime = float(data[0][jj].split()[0])
	
	current_time = 0
	
	for ii in range(1, len(data)):
	    next_file = 1
	    for jj in range(0, len(data[ii])):
		if (len(data[ii][jj].split()) > 0 and data[ii][jj][0] != '#'):
		    current_time = float(data[ii][jj].split()[0])
		    if (current_time > lasttime or next_file == 0):
			fp.write(data[ii][jj])
		else:
		    fp.write(data[ii][jj])
		
		#if (float(data[ii-1][-1].split()[0]) < float(data[ii][jj].split()[0])):
		    #for kk in range(0, len(data[ii][jj])):
			#fp.write(data[ii][jj][kk]+" ")
		    #fp.write("\n")
		#    fp.write(data[ii][jj])
		if (current_time > lasttime):
		    next_file = 0
		    if (len(data[ii][jj].split()) > 0 and data[ii][jj][0] != '#'):
			lasttime = float(data[ii][jj].split()[0])
	
	fp.close()	
