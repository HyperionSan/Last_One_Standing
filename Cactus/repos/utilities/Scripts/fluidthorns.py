#!/usr/bin/env python
# encoding: utf-8
"""
fluidthorns.py

This script interfaces with FluidDB to store and retrieve metadata about
Cactus thorns. Currently requires cclparse.py and Fom (Fluid Object Mapper).

To get Fom run
    
    easy_install Fom 

Created by Eric Seidel on 2010-10-03.
Copyright (c) 2010 . All rights reserved.
"""

__author__ = 'Eric Seidel'

import sys
import os
from pprint import pprint
from optparse import OptionParser
from getpass import getpass
try:
    from fom.dev import sandbox_fluid
    from fom.session import Fluid
    from fom.mapping import Namespace, Object, Tag, tag_value, tag_relation
except ImportError:
    print 'fluidthorns.py requires Fom (Fluid Object Mapper) to run. Please install Fom and try again.'
    sys.exit(1)

# local imports
import cclparse
import parse_readmes

parser = OptionParser()
parser.add_option('-u', '--user', dest='user', default=None, help='login to FluidDB as user USER', metavar='USER')
parser.add_option('-p', '--pass', dest='password', default=None, help='login to FluidDB with password PASS', metavar='PASS')

def main():
    (options, args) = parser.parse_args()
    fluid = None
    user, password = options.user, options.password
    if user:
        fluid = Fluid()
        if password is None:
            password = getpass(prompt="Password for %s: " % user)
    else:
        fluid = sandbox_fluid()
        user = 'test'
        password = 'test'
    fluid.login(user, password)
    fluid.bind()
    ns = Namespace(user)
    try:
        ns.create_namespace(u'CCTK', u'namespace for cactus thorns')
    except:
        # namespace already exists
        pass
    ns = Namespace(u'%s/CCTK' % user)
    
    tags = [u'name', u'description', u'authors', u'language', u'version',
            u'url', u'scm', u'implements', u'provides_function', 
            u'schedules_function', u'requires_function', u'uses_function',
            u'shares', u'inherits', u'arrangement', u'licence',
            u'maintainers']
    
    for tag in tags:
        try:
            ns.create_tag(tag, tag, indexed=False)
            print "%s created" % tag
        except:
            print "%s exists" % tag
    
    thorns = cclparse.parse_deps()
    metas  = parse_readmes.parse()
    for t in thorns:
        print "Creating %s" % t['name']
        thorn = Thorn(about="CCTK:%s" % t['name'])
        # try:
        #     thorn = Thorn(about="CCTK:%s" % t['name'])
        # except:
        #     thorn = Thorn.filter('fluiddb/about="CCTK:%s"' % t['name'])[0]
        #     print "Exists: %s" % t['name']
        thorn.name                  = os.path.basename(t['name'])
        thorn.arrangement           = os.path.dirname(t['name'])
        thorn.implements            = t['implements']
        try:
            thorn.inherits          = t['inherits']
        except KeyError:
            thorn.inherits          = []
        try:
            thorn.provides_function     = t['provides_function']
        except KeyError:
            thorn.provides_function     = []
        #try:
        #    print t['schedules_function']
        #    thorn.schedules_function    = t['schedules_function']
        #except KeyError:
        #    pass
        try:
            thorn.requires_function     = t['requires_function']
        except KeyError:
            thorn.requires_function     = []
        try:
            thorn.uses_function         = t['uses_function']
        except KeyError:
            thorn.uses_function         = []
        try:
            thorn.shares                = t['shares']
        except KeyError:
            thorn.shares                = []
            
        
        # now do meta data
        if t['name'] in metas:
            #print thorn.name
            meta = metas[t['name']]
            #pprint(meta)
            try:
                thorn.authors = meta['authors']
            except KeyError:
                thorn.authors = []
            try:
                thorn.maintainers = meta['maintainers']
            except KeyError:
                thorn.maintainers = []
            try:
                thorn.licence = meta['licence']
            except KeyError:
                thorn.licence = 'Unknown'
            try:
                thorn.description = meta['purpose']
            except KeyError:
                thorn.description = 'Unknown'
            try:
                thorn.scm = meta['scm']
            except KeyError:
                thorn.scm = 'Unknown'
            try:
                thorn.url = meta['url']
            except KeyError:
                thorn.url = 'Unknown'
            try:
                thorn.language = meta['language']
            except KeyError:
                thorn.language = 'Unknown'
            try:
                thorn.version = meta['version']
            except KeyError:
                thorn.version = 'Unknown'
        else:
            # unable to parse README
            thorn.authors = []
            thorn.maintainers = []
            thorn.licence = 'Unknown'
            thorn.description = 'Unknown'
            thorn.scm = 'Unknown'
            thorn.url = 'Unknown'
            thorn.language = 'Unknown'
            thorn.version = 'Unknown'
        thorn.save()
    
    print "The following thorns implement the 'driver' interface:"    
    for o in Thorn.filter('%s/CCTK/implements="driver"' % user):
        print o.about
        print "%s/%s" % (o.arrangement, o.name)
        print "Authors: %s" % o.authors
        print ''

class MetaData(Object):
    """Mixin to provide metadata for a Cactus Thorn."""
    name                = tag_value(u'gridaphobe/CCTK/name')
    arrangement         = tag_value(u'gridaphobe/CCTK/arrangement')
    description         = tag_value(u'gridaphobe/CCTK/description')
    authors             = tag_value(u'gridaphobe/CCTK/authors')
    maintainers         = tag_value(u'gridaphobe/CCTK/maintainers')
    licence             = tag_value(u'gridaphobe/CCTK/licence')
    language            = tag_value(u'gridaphobe/CCTK/language')
    version             = tag_value(u'gridaphobe/CCTK/version')
    url                 = tag_value(u'gridaphobe/CCTK/url')
    scm                 = tag_value(u'gridaphobe/CCTK/scm')
    implements          = tag_value(u'gridaphobe/CCTK/implements')
    provides_function   = tag_value(u'gridaphobe/CCTK/provides_function')
    schedules_function  = tag_value(u'gridaphobe/CCTK/schedules_function')

class Dependency(Object):
    """Mixin to provide depencencies for a Cactus Thorn."""
    requires_function   = tag_value(u'gridaphobe/CCTK/requires_function')
    uses_function       = tag_value(u'gridaphobe/CCTK/uses_function')
    shares              = tag_value(u'gridaphobe/CCTK/shares')
    inherits            = tag_value(u'gridaphobe/CCTK/inherits')

class Thorn(MetaData, Dependency):
    """A Cactus Thorn."""
    def __repr__(self):
        """docstring for __repr__"""
        return 'Thorn: %s' % self.name
    
#####################################################################
if __name__ == '__main__':
    main()
