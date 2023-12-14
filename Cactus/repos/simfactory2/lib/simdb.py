# -*- coding: ascii -*-
# 
# MDB -- MachineDatabase wrapper around pyini.
# Convienient access to the machine database.
# 
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
# 
#

import pyini
import os
import sys

from libutil import *
import simenv,simlib

class Entry:
    def __init__(self, name, dict):
        keys = dict.keys()
        
        self.InternalDictionary = dict
        self.EntryName = name
        
        for key in keys:
            io = dict[key]
            setattr(self, key, io.Value)
    
    def GetKeys(self):
        return self.InternalDictionary.keys()
    
    def HasKey(self, key):
        return hasattr(self, key)
    
    # can also just be accessed as MachineEntry.attribute
    def GetKey(self, key):
        if self.HasKey(key):
            return getattr(self, key)
        
        return None
        
    def HasKeys(self, keys):
        for k in keys:
            if not(self.HasKey(k)):
                return False
        
        return True
        
class ConfigurationDatabase:

    def __init__(self, mdbDirectory=None, cdb=None, udb=None):
        self.mdbDirectory = mdbDirectory
        self.mdbFilename = None
        self.cdbFilename = cdb
        self.udbFilename = udb
        
        self.SyntaxFile = "%s%s%s" % (simenv.SYNTAX_PATH, os.sep, "mdb-syntax.ini")
        self.DescriptionSyntaxFile = "%s%s%s" % (simenv.SYNTAX_PATH, os.sep, "description-syntax.ini")
        self.MachineParser = None
        self.ConfigParser = None
        self.syntaxChecker = None
        self.MachineCache = dict()
        self.ConfigCache = dict()
        
    def Load(self, mdbDirectory=None, cdb=None, udb=None):
        
        if mdbDirectory is not None:
            self.mdbDirectory = mdbDirectory
        
        if cdb is not None:
            self.cdbFilename = cdb
            
        if udb is not None:
            self.udbFilename = udb
            
        if (mdbDirectory is None and self.mdbDirectory is None) or (not(os.path.exists(self.mdbDirectory))):
            fatal("initializing machine database: No database file provided, or is not readable")
        
        if cdb is None and self.cdbFilename is None:
            fatal("cannot initialize configuration database, no database file provided")
        
        # first, verify the correctness of our mdb syntax file.
        sc = pyini.IniSyntaxChecker(self.DescriptionSyntaxFile, self.SyntaxFile)
        sc.SyntaxCheck()
        
        # now load our passed in mdb database
        self.MachineParser = pyini.IniParser()
        
        for ff in os.listdir(self.mdbDirectory):
            
            if not(ff.endswith(".ini")):
                continue
                
            filePath = simlib.BuildPath(simenv.MDB_PATH, ff)

            self.MachineParser.UpdateFromIni(filePath, True)
        
        # load the cdb database, which at the moment has no syntax file and convert any blocks to lists
        self.ConfigParser = pyini.IniParser(self.cdbFilename)
        self.ConvertBlocks()
        
        udbOnlyMachines = []
        if self.udbFilename is not None:
            # import new sections to the machine database, but not the config database.
            nonUdbOnlyMachines = list(self.MachineParser.GetSections())
            self.MachineParser.UpdateFromIni(self.udbFilename, True)
            for m in self.MachineParser.GetSections():
                if m not in nonUdbOnlyMachines:
                    udbOnlyMachines.append(m)
            self.ConfigParser.UpdateFromIni(self.udbFilename)
    
        self.MachineParser.UpdateFromDict(simenv.OptionsManager.MDBKeys)

        # synatax checker
        msc = pyini.IniSyntaxChecker(self.SyntaxFile, self.MachineParser)

        # fill out possibly missing entries from defaults
        msc.SetSectionDefaults()

        # until an udb exists, username and email are (likely) not set so skip checking
        if self.udbFilename is not None:
            for m in list(self.MachineParser.GetSections()):
                # drop all machines that were only in udb but are not complete
                if m in udbOnlyMachines:
                    if not msc.SyntaxCheck(m, lambda t: True):
                        warning("Removing incomplete machine '%s' found only in %s." % (m, self.udbFilename))
                        self.MachineParser.DeleteSection(m)
                else:
                    msc.SyntaxCheck(m)
    
    def ConvertBlocks(self):
        keys = self.ConfigParser.GetGlobalKeys()
        
        for k in keys:
            io = self.ConfigParser.GetGlobalOption(k)
            if io.IsBlock:
                io.ConvertToList()

    # --- ACCESSOR METHODS ---
    
    def GetMachine(self, key):
        cache = self.MachineCache
        if key in cache:
            return cache[key]
        
        parser = self.MachineParser
        if parser.HasSection(key):
            sdict = parser.GetSectionAsDict(key)
            entry = Entry(key, sdict)
            cache[key] = entry
            return entry
        else:
            fatal("Could not find machine entry for '%s'" % (key))
    
    def HasMachineEntry(self, key):
        cache = MachineCache
        if key in cache:
            return True
        
        parser = self.MachineParser
        return parser.HasSection(key)

    def GetMachines(self):
        parser = self.MachineParser
        return parser.GetSections()
        
    def GetConfigOption(self, option):
        
        ret = self.ConfigParser.GetGlobalOption(option)
        
        if ret is not None:
            return self.ConfigParser.GetGlobalOption(option).Value
        else:
            return None
            
    def HasConfigOption(self, option):
        return self.ConfigParser.HasOption(None, option)

    def GetConfigOptions(self):
        return self.ConfigParser.GetGlobalKeys()
