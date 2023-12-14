#!/usr/bin/env python
# encoding: utf-8

'''
parse_readmes.py

This script parses the READMEs of Cactus Arrangements
'''

import os
import re
import sys
from subprocess import *

from pprint import pprint

TITLE_RE        = re.compile(r'^Cactus\s*Code\s*Thorn\s*(\w+)')
AUTHOR_RE       = re.compile(r'Author\(s\)\s*:\s*(.*?)^[^\s]', re.DOTALL | re.MULTILINE)
INDIV_RE        = re.compile(r'^\s*(.+?)\s*$', re.MULTILINE)
MAINTAINER_RE   = re.compile(r'Maintainer\(s\)\s*:\s*(.*?)^[^\s]', re.DOTALL | re.MULTILINE)
LICENCE_RE      = re.compile(r'Licence\s*:\s*(.*)\n')
PURPOSE_RE      = re.compile(r'1.\s*Purpose\n\n(.*?)\n[\n$]', re.DOTALL | re.MULTILINE)

SVN_RE          = re.compile(r'^URL:\s*(.*)$', re.MULTILINE)

def main():
    pprint(parse())

def parse():
    thorns = {}
    for (dirpath, dirnames, filenames) in os.walk('arrangements', followlinks = True):
        if 'README' in filenames and len(dirpath.split(os.sep)) == 3:
            path = os.path.join(dirpath, 'README')
            f = open(path)
            readme = f.read()
            f.close()
        
            thorn = {}
        
            title = TITLE_RE.search(readme)
            if title:
                title = title.group(1)
                thorn['title'] = title
            else:
                print "invalid title in %s" % path
                continue
        
            authors = AUTHOR_RE.search(readme)
            if authors:
                authors = INDIV_RE.findall(authors.group(1))
                thorn['authors'] = authors
            else:
                thorn['authors'] = []
        
            maintainers = MAINTAINER_RE.search(readme)
            if maintainers:
                maintainers = INDIV_RE.findall(maintainers.group(1))
                thorn['maintainers'] = maintainers
            else:
                thorn['maintainers'] = []
            
            licence = LICENCE_RE.search(readme)
            if licence:
                licence = licence.group(1)
                thorn['licence'] = licence
            else:
                thorn['licence'] = None
        
            purpose = PURPOSE_RE.search(readme)
            if purpose:
                purpose = re.sub(r'\n', ' ', purpose.group(1))
                thorn['purpose'] = purpose
            else:
                thorn['purpose'] = None
            
            thorn.update(find_scm_url(dirpath))
            
            canonical_name = os.sep.join(dirpath.split(os.sep)[-2:])
            thorns[canonical_name.lower()] = thorn
    return thorns

def find_scm_url(thorn_path):
    items = os.listdir(thorn_path)
    # determine which scm is being used
    if 'CVS' in items:
        return cvs_url(thorn_path)
    elif '.svn' in items:
        return svn_url(thorn_path)
    elif call(['cd', thorn_path, '&&', 'git', 'branch']) == 0:
        return git_url(thorn_path)
    elif call(['cd', thorn_path, '&&', 'hg', 'branches']) == 0:
        return hg_url(thorn_path)
    else: return {}
        
def cvs_url(thorn_path):
    f = open(thorn_path + os.path.join('CVS', 'Root'))
    url = f.read().strip()
    f.close()
    return {'scm' : 'cvs', 'url' : url}

def svn_url(thorn_path):
    out = Popen(['svn', 'info', thorn_path], stdout=PIPE).communicate()[0]
    url = SVN_RE.search(out).group(1)
    return {'scm' : 'svn', 'url' : url}

def git_url(thorn_path):
    url = Popen(['cd %s && git config --get remote.origin.url' % thorn_path], stdout=PIPE, shell=True).communicate()[0].strip()
    return {'scm' : 'git', 'url' : url}
    
def hg_url(thorn_path):
    url = Popen(['cd %s && hg showconfig paths.default' % thorn_path], stdout=PIPE, shell=True).communicate()[0].strip()
    return {'scm' : 'hg', 'url' : url}

#####################################################
if __name__ == '__main__':
    main()