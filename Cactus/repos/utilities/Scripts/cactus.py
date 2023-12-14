#! /usr/bin/env python
# -*-Python-*-

# Driver for the Cactus software framework
# Copyright (C) 2010 the Cactus Team <cactusmaint@cactuscode.org>

import os
import re
import shutil
import sys
import distutils.dep_util
from optparse import OptionParser
from types import ListType



################################################################################
###   Setup, options, environment   ############################################
################################################################################



# Version number
CCTK_VERSION_MAJOR = '4'
CCTK_VERSION_MINOR = '0'
CCTK_VERSION_OTHER = 'b16'

CCTK_VERSION = CCTK_VERSION_MAJOR + '.' + CCTK_VERSION_MINOR + '.' + CCTK_VERSION_OTHER
os.environ['CCTK_VERSION_MAJOR'] = CCTK_VERSION_MAJOR
os.environ['CCTK_VERSION_MINOR'] = CCTK_VERSION_MINOR
os.environ['CCTK_VERSION_OTHER'] = CCTK_VERSION_OTHER
os.environ['CCTK_VERSION'] = CCTK_VERSION



# Dividers to make the screen output slightly nicer
DIVIDER = '--------------------------------------------------------------------------------'



# Find out where we are
CCTK_HOME = os.getcwd()
os.environ['CCTK_HOME'] = CCTK_HOME

# Find out where the configuration directory is
CONFIGS_DIR = os.environ.get('CACTUS_CONFIGS_DIR', CCTK_HOME + '/configs')
os.environ['CONFIGS_DIR'] = CONFIGS_DIR

# Find all existing configurations
try:
    CONFIGURATIONS = os.listdir(CONFIGS_DIR)
except OSError,e:
    # In case there are no configs yet
    CONFIGURATIONS = []

# Find out how to call the driver
DRIVER = sys.argv[0]



# Parse options
parser = OptionParser()
parser.add_option(''  , '--prompt'   , type='string', dest='prompt'   )
parser.add_option(''  , '--silent'   , type='string', dest='silent'   )
parser.add_option('-j', '--jobs'     , type='int'   , dest='jobs'     )
parser.add_option(''  , '--fjobs'    , type='int'   , dest='fjobs'    )
parser.add_option(''  , '--tjobs'    , type='int'   , dest='tjobs'    )
parser.add_option(''  , '--options'  , type='string', dest='options'  )
parser.add_option(''  , '--thornlist', type='string', dest='thornlist')
parser.add_option(''  , '--thorns'   , type='string', dest='thorns'   )
parser.add_option(''  , '--buildlist', type='string', dest='buildlist')
parser.add_option(''  , '--utils'    , type='string', dest='utils'    )
(options, args) = parser.parse_args()



if options.prompt:
   print 'Option --prompt not supported.'
   exit(1)

if options.silent == None or options.silent.lower() == 'yes':
   SILENT = True
   SILENTMAKE = ['-s']
elif options.silent.lower() == 'no':
   SILENT = False
   SILENTMAKE = []
else:
   print 'Option --silent has wrong value (should be "yes" or "no").'
   exit(1)

if not options.jobs:
   PARFLAGS = []
else:
   PARFLAGS = ['-j', str(options.jobs)]

if not options.fjobs:
   FPARFLAGS = []
else:
   FPARFLAGS = ['-j', str(options.fjobs)]
os.environ['FPARFLAGS'] = ' '.join(FPARFLAGS)

if not options.tjobs:
   TPARFLAGS = []
else:
   TPARFLAGS = ['-j', str(options.tjobs)]

if not options.options:
   SETUP_OPTIONS = []
else:
   SETUP_OPTIONS = ['-config_file=' + options.options]

THORNLIST = options.thornlist

# Directory for configuration options

# Set THORNLIST_DIR to '.' if it is not set already, and if THORNLIST
# does not contain an absolute patch (starting with a slash)
if THORNLIST and THORNLIST[0] == '/':
   THORNLIST_DIR = ''
else:
   THORNLIST_DIR = os.environ.get('THORNLIST_DIR', '.')

if options.thorns:
   print 'Option --thorns not supported.'
   exit(1)

if not options.buildlist:
   BUILDLIST = []
else:
   BUILDLIST = options.buildlist.split()

if not options.utils:
   UTILS = []
else:
   UTILS = options.utils.split()



################################################################################
###   Helper functions   #######################################################
################################################################################

# call a subprocess and capture the output
def pcall(args,output):
    p = os.popen(args,'r')
    for line in p.readlines():
        output += line
    return p.close()

# Convert list items to strings, and quote them so that they are safe
# for the shell
def quotesafe(item):
   return "'" + re.sub("'", "\\'", str(item)) + "'"

# Convert a nested list of arguments into a command suitable for shell
# execution
def makecmd(args):
   # Flatten nested list recursively
   def flatten(item):
      # Remove one level of nesting
      def expand(list):
         return reduce(lambda a,b: a+b, list, [])
      if isinstance(item, ListType):
         return expand(map(flatten,item))
      return [item]
   return ' '.join(map(quotesafe, flatten(args)))



# Output a log message, unless in silent mode
def log(msg):
   if not SILENT:
      print msg

# Execute a command, which is specified by a (possibly nested) list of
# arguments
def execute(args):
   cmd = makecmd(args)
   log ('Executing: ' + cmd)
   return os.system (cmd)

# Delete a file that may or may not exist
def deletefile(file):
   try:
      os.remove(file)
   except OSError:
      pass



################################################################################
###   Commands   ###############################################################
################################################################################



def default():
   # Default action does nothing
   if not CONFIGURATIONS:
      print 'No configurations exist.'
      print '   "' + DRIVER + ' config <name>" will create a configuration called <name>.'
      print '   "' + DRIVER + ' build <name>" will then build this configuration.'
   else:
      print 'The following configurations exist:'
      print '   ' + ', '.join(CONFIGURATIONS) + '.'
      print 'The command "' + DRIVER + ' build <name>" will build a configuration.'
   print 'Use "' + DRIVER + ' help" to see all available commands.'
   print DIVIDER
   exit()



def help(args):
   print '**********************************'
   print '*** Welcome to the Cactus Code ***'
   print '**********************************'
   print
   if not CONFIGURATIONS:
      print 'No configurations exist yet.'
      print '   "' + DRIVER + ' configure <name>" will create a configuration called <name>.'
      print '   "' + DRIVER + ' build <name>" will then build this configuration.'
   else:
      print 'The following configurations exist:'
      print '   ' + ', '.join(CONFIGURATIONS) + '.'
      print 'The command "' + DRIVER + ' build <name>" will build a configuration.'
   print DIVIDER
   print 'There are a range of commands available to act on a configuration.'
   print 'These are called by "' + DRIVER + ' <command> <name>"'
   print 'Valid commands are:'
   print '   build      : build individual thorns of a configuration'
   print '   clean      : clean a configuration'
   print '                (delete all object and dependency files in the configuration)'
   print '   cleandeps  : clean a configuration\'s dependency files'
   print '   cleanobjs  : clean a configuration\'s object files'
   print '   configure  : create a new configuration, or reconfigure an existing one'
   print '                (overwrites previous configuration options)'
   print '   config-info: display the configuration options for a configuration'
   print '   delete     : delete a configuration'
   print '   editthorns : edit the ThornList file'
   print '   examples   : copy thorn parameter files to examples directory'
   print '   realclean  : restore a configuration to an almost new state'
   print '                (delete everything but the config-data directory'
   print '                 and the ThornList file)'
   print '   rebuild    : rebuild a configuration (force the CST to be re-run)'
   print '   reconfigure: reconfigure an existing configuration'
   print '                using its previous configuration options'
   print '   testsuite  : run the test suites'
   print '   ThornGuide : create the thorn manual for a specific configuration'
   print '   thornlist  : overwrite the ThornList file with a list of all existing thorns'
   print '   update     : update the source files for a specific configuration'
   print '   utils      : build a configuration\'s utility programs'
   print '                from repositories'
   print DIVIDER
   print 'There are commands which act on thorns.'
   print 'These are called by "' + DRIVER + ' <command> <thorn>"'
   print 'Valid commands are:'
   print '   ThornDoc: produce the documentation for the thorn'
   print '             in doc/ThornDoc/<arrangement>/<thorn>'
   print DIVIDER
   print 'There are commands to act on arrangements.'
   print 'These are called by "' + DRIVER + ' <command> <arrangement>"'
   print 'Valid commands are:'
   print '   ArrangementDoc: produce documentation for the arrangement'
   print '                   in doc/ArrangementDoc/<arrangement>'
   print '   ThornDocs     : produce the documentation for the thorns'
   print DIVIDER
   print 'There are also the following stand-alone commands:'
   print '   AllDoc         : build all documentation'
   print '   ArrangementDoc : create documentation for all arrangements in'
   print '                    doc/ArrangementDoc'
   print '   checkout       : check out public arrangements/thorns from repositories'
   print '   default        : create a new configuration with a default name'
   print '   diff           : show differences between installed Cactus and the version'
   print '   distclean      : delete all existing configurations'
   print '   downsize       : remove non-essential files in the repositories'
   print '   MaintGuide     : create the maintainers manual doc/MaintGuide.pdf'
   print '   newthorn       : create a new thorn'
   print '   ReferenceManual: create reference manual doc/ReferenceManual.pdf'
   print '   status         : report on the status of the installed Cactus'
   print '                    (when installed from repositories)'
   print '   TAGS           : create an emacs TAGS file'
   print '   tags           : create a vi TAGS file'
   print '   ThornGuide     : create the thorn manual doc/ThornGuide.pdf'
   print '   ThornDoc       : create documentation for all thorns in doc/ThornDoc'
   print '   thorninfo      : give information about all available thorns'
   print '   update         : update flesh and arrangements from repositories'
   print '   UsersGuide     : create the users manual doc/UsersGuide.pdf'
   print DIVIDER

def UsersGuide():
    print DIVIDER
    print "Creating user documentation UsersGuide.pdf"
    os.chdir("doc/UsersGuide")
    print "  Running pdflatex...."
    output = ''
    r=pcall("pdflatex -interaction=nonstopmode UsersGuide.tex",output)
    if r == 0:
        r = pcall("pdflatex -interaction=nonstopmode UsersGuide.tex",output)
    if r == 0:
        r = pcall("pdflatex -interaction=nonstopmode UsersGuide.tex",output)
    err = re.compile("(?m)^!")
    warn = re.compile("(?m)^LaTeX Warning:")
    if err.search(output):
      print "  Problem in UsersGuide.  See doc/UsersGuide/LATEX_MESSAGES."
      exit(1)
    elif warn.search(output):
      print "  For more information see doc/UsersGuide/LATEX_MESSAGES."
    os.rename("UsersGuide.pdf",CCTK_HOME+"/doc/UsersGuide.pdf")
    print "  UsersGuide.pdf created in doc directory."
    print "  Done."
    print DIVIDER


def version(args):
   print 'Cactus ' + CCTK_VERSION
   print DIVIDER



def configure(args):
   if len(args) <= 1:
      print 'No configuration specified.'
      print DIVIDER
      exit(1)
   CONFIGURATION = args[1]
   if not (CONFIGURATION in CONFIGURATIONS):
      print 'Setting up new configuration ' + CONFIGURATION + '.'
   else:
      print 'Reconfiguring existing configuration ' + CONFIGURATION + '.'
      deletefile (CONFIGS_DIR + '/' + CONFIGURATION + '/config-data/make.thornlist')
   if THORNLIST and not os.path.exists(THORNLIST_DIR + '/' + THORNLIST):
      print 'ThornList ' + THORNLIST_DIR + '/' + THORNLIST + ' does not exist.'
      exit(1)
   status = execute (['perl',
                      '-s', 'lib/make/setup_configuration.pl',
                      SETUP_OPTIONS,
                      CONFIGURATION])
   if status != 0:
      print
      print 'Error creating configuration ' + CONFIGURATION + '.'
      deletefile (CONFIGS_DIR + '/' + CONFIGURATION + '/config-data/cctk_Config.h')
      exit(1)
   if THORNLIST:
      shutil.copy(THORNLIST_DIR + '/' + THORNLIST,
                  CONFIGS_DIR + '/' + CONFIGURATION + '/ThornList')
   print DIVIDER
   # print 'Use "' + DRIVER + ' build ' + CONFIGURATION + '" to build the configuration.'
   # print DIVIDER
   build(args)



def reconfigure(args):
   if len(args) <= 1:
      print 'No configuration specified.'
      print DIVIDER
      exit(1)
   CONFIGURATION = args[1]
   if not (CONFIGURATION in CONFIGURATIONS):
      print 'Configuration ' + CONFIGURATION + ' does not exist.'
      print 'Reconfiguration aborted.'
      exit(1)
   if not os.path.exists(CONFIGS_DIR + '/' + CONFIGURATION + '/config-info'):
      print 'Error reconfiguring ' + CONFIGURATION + ': configuration is incomplete.'
      print 'Use "' + DRIVER + ' configure ' + CONFIGURATION + '"to configure the configuration.'
      exit(1)
   status = execute (['perl',
                      '-s', 'lib/make/setup_configuration.pl',
                      SETUP_OPTIONS,
                      CONFIGURATION])
   if status != 0:
      print
      print 'Error reconfiguring ' + CONFIGURATION + '.'
      deletefile (CONFIGS_DIR + '/' + CONFIGURATION + '/config-data/cctk_Config.h')
      exit(1)
   deletefile (CONFIGS_DIR + '/' + CONFIGURATION + '/config-data/make.thornlist')
   print DIVIDER
   # print 'Use "' + DRIVER + ' build ' + CONFIGURATION + '" to build the configuration.'
   # print DIVIDER
   build(args)



def build(args):
   if len(args) <= 1:
      print 'The following configurations exist:'
      print '   ' + ', '.join(CONFIGURATIONS) + '.'
      print 'Please specify a configuration to build.'
      print DIVIDER
      exit(1)
   CONFIGURATION = args[1]
   if not (CONFIGURATION in CONFIGURATIONS):
      print 'Configuration ' + CONFIGURATION + ' does not exist.'
      print 'Build aborted.'
      exit(1)
   if not os.path.exists(CONFIGS_DIR + '/' + CONFIGURATION + '/config-data/cctk_Config.h'):
      print 'Error: Configuration ' + CONFIGURATION + ' is incomplete.'
      print 'Please check the files in "' + CONFIGS_DIR + '/' + CONFIGURATION + '/config-data" for error messages.'
      print 'You can try again to configure using "' + DRIVER + ' configure ' + CONFIGURATION + '"'
      print 'or delete this configuration with "' + DRIVER + ' delete ' + CONFIGURATION + '"'
      print DIVIDER
      exit(1)
   if distutils.dep_util.newer (CCTK_HOME + '/lib/make/force-reconfigure',
                                CONFIGS_DIR + '/' + CONFIGURATION + '/config-info'):
      print 'Error: Configuration ' + CONFIGURATION + ' is out of date.'
      print '   Please reconfigure your configuration by running the command'
      print
      print '      ' + DRIVER + ' reconfigure ' + CONFIGURATION
      print
      print '   (It is likely that recent changes to the flesh require this.)'
      print DIVIDER
      exit(1)
   execute (['make',
             SILENTMAKE,
             '-f', CCTK_HOME + '/lib/make/make.configuration',
             'TOP=' + CONFIGS_DIR + '/' + CONFIGURATION,
             'CCTK_HOME=' + CCTK_HOME,
             PARFLAGS, TPARFLAGS,
             'rebuild'])
   if BUILDLIST:
      print DIVIDER
      print "Building thorns [" + ', '.join(BUILDLIST) + "] of configuration " + CONFIGURATION + ':'
      print DIVIDER
      os.environ['BUILDLIST'] = ' '.join(BUILDLIST)
      # os.chdir (CONFIGS_DIR + '/' + CONFIGURATION)
      execute (['make',
                SILENTMAKE,
                '-f', CCTK_HOME + '/lib/make/make.configuration',
                'TOP=' + CONFIGS_DIR + '/' + CONFIGURATION,
                'CCTK_HOME=' + CCTK_HOME,
                'build'])
   else:
      print DIVIDER
      print 'Building configuration ' + CONFIGURATION + ':'
      print DIVIDER
      execute (['make',
                SILENTMAKE,
                '-f', CCTK_HOME + '/lib/make/make.configuration',
                'TOP=' + CONFIGS_DIR + '/' + CONFIGURATION,
                'CCTK_HOME=' + CCTK_HOME,
                PARFLAGS, TPARFLAGS])
   print DIVIDER



def rebuild(args):
   if len(args) <= 1:
      print 'The following configurations exist:'
      print '   ' + ', '.join(CONFIGURATIONS) + '.'
      print 'Please specify a configuration to rebuild.'
      print DIVIDER
      exit(1)
   CONFIGURATION = args[1]
   if not (CONFIGURATION in CONFIGURATIONS):
      print 'Configuration ' + CONFIGURATION + ' does not exist.'
      print 'Rebuild aborted.'
      exit(1)
   deletefile (CONFIGS_DIR + '/' + CONFIGURATION + '/config-data/make.thornlist')
   build(args)



def utils(args):
   if len(args) <= 1:
      print 'The following configurations exist:'
      print '   ' + ', '.join(CONFIGURATIONS) + '.'
      print 'Please specify a configuration for which the utilities should be built.'
      print DIVIDER
      exit(1)
   CONFIGURATION = args[1]
   if not (CONFIGURATION in CONFIGURATIONS):
      print 'Configuration ' + CONFIGURATION + ' does not exist.'
      print 'Building utilities aborted.'
      exit(1)
   print 'Building utilities for ' + CONFIGURATION + ':'
   # os.chdir (CONFIGS_DIR + '/' + CONFIGURATION)
   execute (['make',
             SILENTMAKE,
             '-f', CCTK_HOME + '/lib/make/make.configuration',
             'TOP=' + CONFIGS_DIR + '/' + CONFIGURATION,
             'CCTK_HOME=' + CCTK_HOME,
             PARFLAGS, TPARFLAGS,
             'utils',
             'UTILS=' + ' '.join(UTILS),
             'CONFIG_NAME=' + CONFIGURATION])
   print DIVIDER



def config_info(args):
   if len(args) <= 1:
      print 'The following configurations exist:'
      print '   ' + ', '.join(CONFIGURATIONS) + '.'
      print 'Please specify a configuration for which to show the configuration information.'
      print DIVIDER
      exit(1)
   CONFIGURATION = args[1]
   if not (CONFIGURATION in CONFIGURATIONS):
      print 'Configuration ' + CONFIGURATION + ' does not exist.'
      print 'Displaying configuration information aborted.'
      exit(1)
   print 'Displaying configuration information for ' + CONFIGURATION + ':'
   print DIVIDER
   # TODO: Should use Python mechanism for this
   execute (['cat', 'configs/' + CONFIGURATION + '/config-info'])
   print DIVIDER



#elif COMMAND == 'cleandeps':
#elif COMMAND == 'cleanobjs':
#elif COMMAND == 'clean':
#elif COMMAND == 'realclean':
#
#elif COMMAND == 'delete':




def error(args):
   print 'Unknown command "' + COMMAND + '".'
   print 'Use "' + DRIVER + ' help" to see all available commands.'
   print DIVIDER
   exit(1)



################################################################################
###   Main   ###################################################################
################################################################################



print DIVIDER
print 'Cactus ' + CCTK_VERSION
print DIVIDER
if not args:
   default()
else:
   COMMAND = args[0]
   if COMMAND == 'help': help(args)
   elif COMMAND == 'UsersGuide': UsersGuide()
   elif COMMAND == 'version': version(args)
   elif COMMAND == 'configure': configure(args)
   elif COMMAND == 'reconfigure': reconfigure(args)
   elif COMMAND == 'build': build(args)
   elif COMMAND == 'rebuild': rebuild(args)
   elif COMMAND == 'utils': utils(args)
   elif COMMAND == 'config-info': config_info(args)
   else: error(args)
