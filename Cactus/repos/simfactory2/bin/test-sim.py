#!/usr/bin/env python3
# -*- coding: ascii -*-
#
# Test the Simulation Factory
# 
# This file has a number of tests
# that are defined as functions.
# Each function name should start
# with the string "test_" and
# should require no arguments.
#
# If a test has run, its output 
# will be in tests/test_name.txt.
# If the test succeeds, it will
# end with the string "PASS".
# If a test has passed, it will
# not be re-run.

from traceback import print_exc
from subprocess import Popen, PIPE
from time import sleep
import sys
try:
    if sys.stdout.isatty():
        from termcolor import colored
except:
    pass

if "colored" not in globals():
    def colored(x,_):
        return x

import os
import re

class SkipTest(Exception):
    def __init_(self, msg):
        Exception.__init__(self, msg)

def get_remote():
    host = os.environ.get("REMOTE", None)
    if host is None:
        raise SkipTest("You need to set the REMOTE variable for tests of sync, etc.")
    return host

log_fd = None

thorn_list = "./simfactory/etc/thornlists/expressions.th"
alt_thorn_list = "alt.th"

def cmd(args,inp="",silent=False,ret_codes=[0]):
    if not silent:
        print("cmd:",args)
        print("cmd:",args,file=log_fd)
    p = Popen(args, universal_newlines=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate(inp)
    if not silent:
        if out.strip() != '':
            print(out, file=log_fd)
        if err.strip() != '':
            print(err, file=log_fd)
    assert p.returncode in ret_codes
    return out, err

# Make sure we have an MPI installed
cmd(["mpirun","--help"],silent=True)

defs = "./simfactory/etc/defs.local.ini"

def test_config1():
    if os.path.exists(defs):
        os.remove(defs)
    out, err = cmd(["./simfactory/bin/sim","setup-silent"])
    assert "Contents successfully written to" in out
    assert os.path.exists(defs)
    g = re.search(r'Determining local machine name: ([\S+]+)', out)
    mach_file = os.path.join("simfactory","mdb", "machines", g.group(1) + ".ini")
    assert os.path.exists(mach_file), mach_file

def test_config1b():
    if os.path.exists(defs):
        os.remove(defs)
    out, err = cmd(["./simfactory/bin/sim","setup-silent","--setup-user",os.getlogin()])
    assert "Contents successfully written to" in out
    assert os.path.exists(defs)
    g = re.search(r'Determining local machine name: ([\S+]+)', out)
    mach_file = os.path.join("simfactory","mdb", "machines", g.group(1) + ".ini")
    assert os.path.exists(mach_file), mach_file

def test_config1c():
    if os.path.exists(defs):
        os.remove(defs)
    out, err = cmd(["./simfactory/bin/sim","setup-silent","--setup-email","foo@bar.com"])
    assert "Contents successfully written to" in out
    assert os.path.exists(defs)
    g = re.search(r'Determining local machine name: ([\S+]+)', out)
    mach_file = os.path.join("simfactory","mdb", "machines", g.group(1) + ".ini")
    assert os.path.exists(mach_file), mach_file

def test_config2():
    if os.path.exists(defs):
        os.remove(defs)
    out, err = cmd(["./simfactory/bin/sim","setup"],inp="\n"*7)
    assert "Contents successfully written to" in out
    assert err.strip() == ""
    assert os.path.exists(defs)
    g = re.search(r'Determining local machine name: ([\S+]+)', out)
    mach_file = os.path.join("simfactory","mdb", "machines", g.group(1) + ".ini")
    assert os.path.exists(mach_file), mach_file

def build_ok(out,exe,sim="sim"):
    assert os.path.exists(exe)
    assert ("Done creating cactus_%s." % sim) in out
    assert "All done !" in out
    assert ("Building utilities for %s" % sim) in out
    assert "Done." in out

exe = os.path.join("exe", "cactus_sim")
def test_build1():
    if os.path.exists(exe):
        os.remove(exe)
    out, err = cmd(["./simfactory/bin/sim","build","--thornlist",thorn_list])
    build_ok(out, exe)

def test_buildj():
    # Test -j option
    cmd(["rm","-fr","configs"])
    out, err = cmd(["./simfactory/bin/sim","build","-j6","--thornlist",thorn_list])
    build_ok(out, exe)

def test_build_mdbkey():
    # Test mdbkey
    cmd(["rm","-f",exe])
    out, err = cmd(["./simfactory/bin/sim","build","--mdbkey","make","make -j6","--thornlist",thorn_list])
    build_ok(out, exe)

def test_build_verbose():
    # Test verbose
    alt_exe = os.path.join("exe", "cactus_alt")
    cmd(["rm","-f",alt_exe])
    out, err = cmd(["./simfactory/bin/sim","build","alt","--verbose","--thornlist",thorn_list])
    build_ok(out, alt_exe, "alt")
    assert "Info: EXECUTING COMMAND" in out

# ppn             = 48
# max-num-threads = 48
# num-threads     = 24
# nodes           = 1
def get_num(name, data):
    g = re.search(name+r'\s*=\s*(\d+)', data)
    return int(g.group(1))

def test_whoami():
    global mach_file, ppn, nodes, home
    home = os.environ["HOME"]
    out, err = cmd(["./simfactory/bin/sim","whoami"],silent=True)
    g = re.search(r"Current machine: (\S+)", out)
    mach_file = os.path.join("simfactory","mdb", "machines", g.group(1) + ".ini")
    with open(mach_file, "r") as fd:
        ini_data = fd.read()
        ppn = get_num("ppn", ini_data)
        nodes = get_num("nodes", ini_data)
    return True

def test_createrun():
    cmd(["rm","-fr",os.path.join(home,"simulations/mytests")])
    out, err = cmd(["./simfactory/bin/sim","create-run","mytests","--testsuite", \
        "--procs",str(ppn),"--ppn-used",str(ppn)])
    assert "Skeleton Created" in out
    assert re.search(r'Number failed\s+->\s+0\b', out)

def test_createrunpar():
    cmd(["rm","-fr",os.path.join(home,"simulations/mytestsP")])
    out, err = cmd(["./simfactory/bin/sim","create-run","mytestsP","--parfile", \
        "repos/cactustest/TestPar/test/expressions.par", \
        "--procs",str(ppn),"--ppn-used",str(ppn)])
    assert ": String 9:" in out

def test_failcreaterun():
    out, err = cmd(["./simfactory/bin/sim","create-run","mytests","--testsuite", \
        "--procs",str(ppn),"--ppn-used",str(ppn)],ret_codes=[1])
    assert "already exists" in err 

def test_create():
    home = os.environ["HOME"]
    cmd(["rm","-fr",os.path.join(home,"simulations/mytests2")])
    out, err = cmd(["./simfactory/bin/sim","create","mytests2","--testsuite"])
    assert "Skeleton Created" in out

def test_run():
    out, err = cmd(["./simfactory/bin/sim","run","mytests2", \
        "--procs",str(ppn),"--ppn-used",str(ppn)])
    assert re.search(r'Number failed\s+->\s+0\b', out)

def test_list_sim():
    out, err = cmd(["./simfactory/bin/sim","list-sim"])
    assert re.search(r'mytests.*ACTIVE.*FINISHED', out)
    assert re.search(r'mytests2.*ACTIVE.*FINISHED', out)

def test_list_conf():
    out, err = cmd(["./simfactory/bin/sim","list-conf"])
    assert re.search(r'\bsim.*built', out)
    assert re.search(r'\balt.*built', out)

def test_create2():
    home = os.environ["HOME"]
    cmd(["rm","-fr",os.path.join(home,"simulations/mytests3")])
    out, err = cmd(["./simfactory/bin/sim","create","mytests3","--testsuite"])
    assert "Skeleton Created" in out

def test_submit():
    out, err = cmd(["./simfactory/bin/sim","submit","mytests3", \
        "--procs",str(ppn),"--ppn-used",str(ppn)])
    sleep(5)

def test_follow():
    out, err = cmd(["./simfactory/bin/sim","show-output","--follow","mytests3"])
    assert re.search(r'Number failed\s+->\s+0\b', out)

def test_checkout():
    with open(thorn_list,"r") as fd:
        contents = fd.read()
        contents = re.sub(r'## CactusTest/TestMath','CactusTest/TestMath',contents)
    with open(alt_thorn_list, "w") as fd:
        fd.write(contents)
    if not os.path.exists("GetComponents"):
        cmd(["curl","-kLO","https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents"])
    os.chmod("GetComponents",0o0755)
    out, err = cmd(["./simfactory/bin/sim","checkout",alt_thorn_list],inp="no\n")
    assert "components checked out successfully" in out
    assert "Time Elapsed:" in out

def test_reconfig():
    # Test verbose
    alt_exe = os.path.join("exe", "cactus_alt")
    out, err = cmd(["./simfactory/bin/sim","build","alt","--reconfig","--thornlist",alt_thorn_list])
    build_ok(out, alt_exe, "alt")
    out, err = cmd(["./exe/cactus_alt","-T"],ret_codes=[0,1])
    assert "TestMath" in out

def test_get_outputdir():
    out, err = cmd(["./simfactory/bin/sim","get-output-dir","mytests"])
    assert "output-0000" in out

def test_purge():
    out, err = cmd(["./simfactory/bin/sim","list-sim"])
    g = re.search(r'(mytests\w+)[ \t]+\[', out)
    out, err = cmd(["./simfactory/bin/sim","purge",g.group(1)])
    assert "has been moved to trash folder" in out
    out, err = cmd(["./simfactory/bin/sim","list-sim"])
    assert not re.search(r'\b'+g.group(1)+r'\b', out)

def test_sync():
    out, err = cmd(["./simfactory/bin/sim","sync","--mdbkey","sshopts","-oStrictHostKeyChecking=no",get_remote()])
    assert "sending incremental file list" in out
    assert re.search(r'total size is.*speedup is', out)

os.makedirs("tests", exist_ok=True)
all_vars = [k for k in globals()]
success = True
for k in all_vars:
    v = globals()[k]
    if re.match(r'^test_.*',k) and type(v) == type(test_list_conf):
        test_path = os.path.join("tests", k+".txt")
        if os.path.exists(test_path):
            with open(test_path,"r") as check_fd:
                passed = check_fd.read().strip().endswith("PASS")
        else:
            passed = False
        if passed:
            print(colored("Already run:","green"),k)
        else:
            print(colored("Running:","yellow"),k,"...")
            log_fd = open(test_path, "w")
            try:    
                ret = v()
                if ret is None:
                    print("PASS", file=log_fd)
                print("  ",colored("** Test succeeded **","green"))
            except SkipTest as s:
                print("  ",colored("** Skipping test **","yellow"))
                print(str(s))
            except:
                print("  ",colored("** Test failed **","red"))
                success = False
                print_exc()
                break
            finally:
                log_fd.close()

if success:
    print(colored("All Tests Passed", "green"))
