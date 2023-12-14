#!/usr/bin/env python

# this script queries an xml list of cactus thorn dependencies,
# and generates a list of required thorns based on the specified thorn

import os
import sys
from pprint import pprint
import re
from xml.dom.minidom import parse
from optparse import OptionParser

REQUIRED_THORNS = {}
dom = parse('cactusdeps.xml')
thorns = dom.getElementsByTagName('Thorn')

usage = "usage: %prog [options] thorn"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--all-functions", action="store_true",
                    dest="all_functions", default=False,
                    help="treated used functions as required functions")
parser.add_option("-d", "--debug", action="store_true",
                    dest="debug", default=False,
                    help="print out debug messages")
parser.add_option("-a", "--all-thorns", action="store_true",
                    dest="all_thorns", default=False,
                    help="find dependencies for all thorns instead of single query")
parser.add_option('-b', '--bin', action='store', dest='timebin',
                    default=None, help='print list of thorns that schedule functions in the specified bin')
(options, args) = parser.parse_args()

def main():
    if options.timebin is not None:
        search_timebins(options.timebin)
        sys.exit(0)

    if options.all_functions is True:
        print "finding all functions!!"

    if options.all_thorns is True:
        for th in thorns:
            find_deps(th.getAttribute('Name'))
    else:
        query = args[0]
        find_deps(query)
    #REQUIRED_THORNS.sort()
    print "=========================================================="
    #print "%s requires:" % query
    for th, deps in REQUIRED_THORNS.iteritems():
        for dep in deps:
            print "%s is required by %s" % (th, dep)

    try:
        fptr = open('out.dot', 'w')
    except IOError:
        print "Could not open out.dot"
        sys.exit(1)

    fptr.write('digraph test {\noverlap=false\nsplines=true\n')
    for th, deps in REQUIRED_THORNS.iteritems():
        for dep in deps:
            fptr.write("%s -> %s\n" % (os.path.basename(dep), os.path.basename(th)))
    fptr.write('}')

    fptr.close()

    sys.exit(0)


def find_deps(thorn):
    if options.debug is True:
        print "Searching for thorn: %s" % thorn
    for th in thorns:
        # search for thorn matching input thorn
        if thorn.lower() == th.getAttribute('Name').lower():
            # loop through inherited thorns
            deps = th.getElementsByTagName('Inherits')
            for dep in deps:
                print "%s inherits from %s" % (thorn, dep.firstChild.data.strip())
                find_thorns(dep.firstChild.data.strip(), thorn)

            # loop through required functions
            reqs = th.getElementsByTagName('Requires_Function')
            for req in reqs:
                func = req.firstChild.data.strip()
                print "%s requires function %s" % (thorn, func)
                find_provides(func, thorn)

            if options.all_functions is True:
                # loop through used functions
                reqs = th.getElementsByTagName('Uses_Function')
                for req in reqs:
                    func = req.firstChild.data.strip()
                    print "%s uses function %s" % (thorn, func)
                    find_provides(func, thorn)


            # loop through required thorns
            reqs = th.getElementsByTagName('Requires_Thorn')
            for req in reqs:
                name = req.firstChild.data.strip()
                print "%s requires thorn %s" % (thorn, name)
                find_req_thorn(name, thorn)

            # loop through thorns sharing variables
            shares = th.getElementsByTagName('Shares')
            for share in shares:
                impl = share.firstChild.data.strip()
                print "%s shares variables with implementation: %s" % (thorn, impl)
                find_thorns(impl, thorn)


def find_thorns(impl, req):
    if options.debug is True:
        print "Searching for implementation: %s" % impl
    # search for implementation matching dep
    implementations = []
    for th in thorns:
        name = th.getAttribute('Name')
        i = th.getAttribute('Implements')
        if impl.lower() == i.lower():
            implementations.append(name)

    if options.all_thorns is True:
        for name in implementations:
            #if not name in REQUIRED_THORNS.keys():
                #REQUIRED_THORNS[name] = req
                # recursive step
                #find_deps(name)
            try:
                if not req in REQUIRED_THORNS[name]:
                    REQUIRED_THORNS[name].append(req)
            except KeyError:
                REQUIRED_THORNS[name] = [req]
                find_deps(name)

    elif len(implementations) > 1:
    	# no need to prompt twice about the same implementation
    	for p in implementations:
    		if p in REQUIRED_THORNS:
    			return

        print "%s is required, but %d thorns implement it." % (impl, len(implementations))
        i = 1
        for imp in implementations:
            print "[%d]  %s" % (i, imp)
            i += 1
        answer = raw_input("--> ")
        name = implementations[int(answer)-1]

        #if not name in REQUIRED_THORNS.keys():
            #REQUIRED_THORNS[name] = req
            # recursive step
            #find_deps(name)
        try:
            if not req in REQUIRED_THORNS[name]:
                REQUIRED_THORNS[name].append(req)
        except KeyError:
            REQUIRED_THORNS[name] = [req]
            find_deps(name)


    elif len(implementations) == 1:
        name = implementations[0]
        #if not name in REQUIRED_THORNS.keys():
            #REQUIRED_THORNS[name] = req
            # recursive step
            #find_deps(name)
        try:
            if not req in REQUIRED_THORNS[name]:
                REQUIRED_THORNS[name].append(req)
        except KeyError:
            REQUIRED_THORNS[name] = [req]
            find_deps(name)
    else:
        pass


def find_req_thorn(thorn, req):
    if options.debug is True:
        print "Searching for thorn: %s" % thorn
    # search for thorn matching thorn
    for th in thorns:
        name = th.getAttribute('Name')
        if re.search(thorn.lower(), name.lower()):
            #if not name in REQUIRED_THORNS.keys():
                #REQUIRED_THORNS[name] = req
                # recursive step
                #find_deps(name)
            try:
                if not req in REQUIRED_THORNS[name]:
                    REQUIRED_THORNS[name].append(req)
            except KeyError:
                REQUIRED_THORNS[name] = [req]
                find_deps(name)

def find_provides(func, req):
    if options.debug is True:
        print "Searching for function: %s" % func
    # search for thorns providing function func
    for th in thorns:
        name = th.getAttribute('Name')
        provides = th.getElementsByTagName('Provides_Function')
        # loop through provided functions
        for prov in provides:
            if func == prov.firstChild.data.strip():
                print "%s provides function %s" % (name, func)
                #if not name in REQUIRED_THORNS.keys():
                    #REQUIRED_THORNS[name] = req
                    # recursive step
                    #find_deps(name)
                try:
                    if not req in REQUIRED_THORNS[name]:
                        REQUIRED_THORNS[name].append(req)
                except KeyError:
                    REQUIRED_THORNS[name] = [req]
                    find_deps(name)

def search_timebins(bin):
    thorns_in_bin = []
    for th in thorns:
        name = th.getAttribute('Name')
        schedules = th.getElementsByTagName('Schedules_Function')
        for s in schedules:
            if len(s.getElementsByTagName(bin)) != 0:
                thorns_in_bin.append(str(name))
                break
    pprint(thorns_in_bin)

###############################################################

if __name__ == '__main__':
    main()