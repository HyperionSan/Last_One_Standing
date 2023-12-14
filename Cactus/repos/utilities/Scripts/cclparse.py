#!/usr/bin/env python
# encoding: utf-8

# this script parses all interface.ccl files located recursively in the current directory
# and generates a database of dependencies

__author__ = 'Eric Seidel'

import os
import sys
from pprint import pprint
import re
from xml.dom.minidom import getDOMImplementation
from optparse import OptionParser

try:
    import yaml
    YAML_FOUND = True
except ImportError:
    #print "Yaml package for python not found on this system. Will only output XML."
    YAML_FOUND = False

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--schedule-chart", action="store_true",
                        dest="schedule_chart", default=False,
                        help="output chart of schedule bins and thorns that schedule functions")
    (options, args) = parser.parse_args()
    
    Dependencies = parse_deps()
    
    if YAML_FOUND is True:
        print yaml.dump(Dependencies)

    to_xml(Dependencies)
    
    if options.schedule_chart is True:
        pprint(schedule_dict)

def parse_deps():
    """Parse the Cactus CCL files and generate dictionary of dependencies"""
    Files = []
    Dependencies = []
    TimeBins = ['recover_parameters', 'startup', 'wragh', 'paramcheck',
                'preregridinitial', 'postregridinitial', 'basegrid', 
                'initial', 'postinitial', 'postrestrictinitial', 
                'postpostinitial', 'recover_variables', 
                'post_recover_variables', 'cpinitial', 'checkpoint', 
                'preregrid', 'postregrid', 'prestep', 'evol', 'postrestrict', 
                'poststep', 'analysis', 'terminate', 'shutdown']

    implement_re = re.compile('implements:\s*(\w+)', re.I)
    inherit_re = re.compile('inherits:\s*(.+)', re.I)
    provides_function_re = re.compile('PROVIDES\s+FUNCTION\s+(\w+)', re.I)
    uses_function_re = re.compile('USES\s+FUNCTION\s+(\w+)', re.I)
    requires_function_re = re.compile('REQUIRES\s+FUNCTION\s+(\w+)', re.I)
    shares_re = re.compile('shares:\s*(\w+)', re.I)
    requires_thorn_re = re.compile('REQUIRES\s+(?!FUNCTION\s*)(\w+)', re.I)
    schedules_function_re = re.compile('schedule\s+(?:group\s+)?(\w+)\s+(?:in|at)\s+(\w+)', re.I)

    # find all interface.ccl and param.ccl files in cwd
    for dirpath, dirnames, filenames in os.walk('arrangements', followlinks=True):
        for file in filenames:
            if file == 'interface.ccl':
                Files.append(os.path.join(dirpath, file))


    for file in Files:
        # first parse interface.ccl
        try:
            fptr = open(file, 'r')
        except IOError:
            print "Could not open %s" % file

        lines = fptr.readlines()

        try:
            fptr.close()
        except IOError:
            print "Could not close %s" % file

        # then parse param.ccl
        file = re.sub('interface.ccl', 'param.ccl', file)

        try:
            fptr = open(file, 'r')
        except IOError:
            print "Could not open %s" % file

        lines += fptr.readlines()

        try:
            fptr.close()
        except IOError:
            print "Could not close %s" % file

        # then configuration.ccl
        file = re.sub('param.ccl', 'configuration.ccl', file)

        try:
            fptr = open(file, 'r')
            lines += fptr.readlines()
            fptr.close()
        except IOError:
            pass

        # then schedule.ccl
        file = re.sub('configuration.ccl', 'schedule.ccl', file)

        try:
            fptr = open(file, 'r')
            lines += fptr.readlines()
            fptr.close()
        except IOError:
            pass

        # get the thorn dir and its parent
        thornname = os.path.basename(os.path.dirname(file))
        parentdir = os.path.basename(os.path.dirname(os.path.dirname(file)))
        thornname = os.path.join(parentdir, thornname)
        file_dict = {'name' : thornname.lower()}
        for line in lines:
            line = line.strip()
            m = re.match(implement_re, line)
            if m:
                file_dict['implements'] = m.group(1).lower()

            m = re.match(inherit_re, line)
            if m:
                inheritance = re.split('\W+', m.group(1).lower())
                file_dict['inherits'] = inheritance

            m = re.match(provides_function_re, line)
            if m:
                try:
                    file_dict['provides_function'].append(m.group(1).lower())
                except KeyError:
                    file_dict['provides_function'] = [m.group(1).lower()]

            m = re.match(uses_function_re, line)
            if m:
                try:
                    file_dict['uses_function'].append(m.group(1).lower())
                except KeyError:
                    file_dict['uses_function'] = [m.group(1).lower()]

            m = re.match(requires_function_re, line)
            if m:
                try:
                    file_dict['requires_function'].append(m.group(1).lower())
                except KeyError:
                    file_dict['requires_function'] = [m.group(1).lower()]

            m = re.match(requires_thorn_re, line)
            if m:
                requires = re.split('\W+', m.group(1).lower())
                # sometimes we have 'REQUIRES THORNS' instead of 'REQUIRES'
                if requires[0].lower() == 'thorns':
                    del requires[0]
                file_dict['requires_thorn'] = requires

            m = re.match(shares_re, line)
            if m:
                try:
                    file_dict['shares'].append(m.group(1).lower())
                except KeyError:
                    file_dict['shares'] = [m.group(1).lower()]

            m = re.match(schedules_function_re, line)
            if m:
                bin, func = m.group(2).lower(), m.group(1).lower()
                if bin in TimeBins:
                    bin = 'cctk_' + bin
                func_dict = {bin : func}
                try:
                    file_dict['schedules_function'].append(func_dict)
                except KeyError:
                    file_dict['schedules_function'] = [func_dict]


        Dependencies.append(file_dict)

    return Dependencies

def to_xml(Dependencies, outfile='cactusdeps.xml'):
    """Write dependencies to an XML file"""
    impl = getDOMImplementation()

    newdoc = impl.createDocument(None, "Cactus", None)
    root = newdoc.documentElement

    schedule_dict = {}

    for dep in Dependencies:
        thorn = newdoc.createElement("Thorn")

        thorn.setAttribute('Name', dep['name'])

        #implements = newdoc.createElement("Implements")
        try:
            thorn.setAttribute('Implements', dep['implements'])
        except KeyError:
            pass

        try:
            for i in dep['inherits']:
                inherits = newdoc.createElement("Inherits")
                s = newdoc.createTextNode(i)
                inherits.appendChild(s)
                thorn.appendChild(inherits)
        except KeyError:
            pass

        try:
            for i in dep['provides_function']:
                provides = newdoc.createElement("Provides_Function")
                s = newdoc.createTextNode(i)
                provides.appendChild(s)
                thorn.appendChild(provides)
        except KeyError:
            pass

        try:
            for i in dep['uses_function']:
                uses = newdoc.createElement("Uses_Function")
                s = newdoc.createTextNode(i)
                uses.appendChild(s)
                thorn.appendChild(uses)
        except KeyError:
            pass

        try:
            for i in dep['requires_function']:
                requires = newdoc.createElement("Requires_Function")
                s = newdoc.createTextNode(i)
                requires.appendChild(s)
                thorn.appendChild(requires)
        except KeyError:
            pass

        try:
            for i in dep['requires_thorn']:
                requires = newdoc.createElement("Requires_Thorn")
                s = newdoc.createTextNode(i)
                requires.appendChild(s)
                thorn.appendChild(requires)
        except KeyError:
            pass

        try:
            for i in dep['shares']:
                shares = newdoc.createElement("Shares")
                s = newdoc.createTextNode(i)
                shares.appendChild(s)
                thorn.appendChild(shares)
        except KeyError:
            pass

        try:
            for i in dep['schedules_function']:
                bin, func = i.keys()[0], i.values()[0]
                schedules = newdoc.createElement("Schedules_Function")
                b = newdoc.createElement(bin)
                s = newdoc.createTextNode(func)
                b.appendChild(s)
                schedules.appendChild(b)
                thorn.appendChild(schedules)

                if options.schedule_chart is True:
                    try:
                        if not func in schedule_dict[bin][dep['name']]:
                            schedule_dict[bin][dep['name']].append(func)
                    except KeyError:
                        try:
                            schedule_dict[bin][dep['name']] = [func]
                        except KeyError:
                            schedule_dict[bin] = {dep['name'] : [func]}
        except KeyError:
            pass

        root.appendChild(thorn)

    try:
        fptr = open(outfile, 'w')
    except IOError:
        print "Could not open %s" % outfile
        sys.exit(1)

    fptr.write(root.toprettyxml(encoding='utf-8'))

    try:
        fptr.close()
    except IOError:
        print "Could not close %s" % outfile
        sys.exit(1)
######################################################################
if __name__ == '__main__':
    main()