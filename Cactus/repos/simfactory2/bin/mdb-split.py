#!/usr/bin/python
# -*- coding: ascii -*-
import os, re

#mdb split

fptr = open('mdb.ini', 'r')

lines = fptr.readlines()

rr = '^\[([0-9a-zA-Z]+)\]$'
rc = re.compile(rr)

machine = None

for line in lines:
	matches = rc.match(line)
	
	if matches is not None:
		if machine is not None:
			print "writing contents to mdb/%s.ini" % machine
			ff = open("mdb/%s.ini" % machine, 'w')
			ff.write(contents)
			ff.close()

		machine = matches.group(1)
		contents = "%s" % line
	else:
		contents = "%s%s" % (contents, line) 

if machine is not None:
	print "writing contents to mdb/%s.ini" % machine
	ff = open("mdb/%s.ini" % machine, 'w')
	ff.write(contents)
	ff.close()
